# redoing it without mediators as covariates
# using minimum distance prison as a proxy for next of kin home

rm(list = ls())
library(SuperLearner)
library(matlabr)
library(MatchIt)
dat <- read.csv('~jacquelinemauro/MergedData.csv')[,-1]
dat <- dat[-which(dat$NCRecid.Event == 'reincarceration'),] #drop if they go back for parole violation

nm = names(dat)
pris.dummies <- dat[,which(names(dat)=='alb'):which(names(dat)=='pit')]
avals <- sort(names(pris.dummies))
pris.dummies = pris.dummies[,avals] #sorts pris.dummies alphabetically so it matches w A
dat$A <- apply(pris.dummies,1,which.max)
dat$A <- factor(avals[dat$A], levels = avals)
dist.df <- dat[,which(nm=='ALB_time'):which(nm=='WAM_time')]
min.prison <- apply(dist.df,1,which.min)
min.prison <- factor(avals[min.prison], levels = avals)

# leave visits and distance out
covs = cbind(dat[,which(nm=="loslastloc"):which(nm=='white')],dat[,which(nm=='urban'):which(nm=='ageyrs')],
             dat[,which(nm=='custody_level'):which(nm=='numofpriorinc')],
             dat[,which(nm=='mh'):which(nm=='highschoolgrad')],
             dat[,which(nm=='numofpriormisconducts')],min.prison
)
# standardize <- function(x){(x-mean(x, na.rm = T))/sd(x, na.rm = T)}
# covsS = apply(covs[,-c(2,3,5,6,11,12,14)],2,standardize)

df <- data.frame(Y = dat$NCRecid3, A = dat$A, covs)
df[,dim(df)[2]] <- factor(df[,dim(df)[2]], levels = sort(unique(df[,dim(df)[2]] )))
to.keep <- complete.cases(df)
df <- df[complete.cases(df),]  # highschool grad the most missing, 63 unobserved values
names(df) <- c('y', 'A',sapply(c(1:dim(covs)[2]), function(k) paste('x',k,sep = "")))
pris.dummies <- pris.dummies[to.keep,]

nms <- c('Recidivism','Prison','Length of Stay', 'White',
         'Urban',"Prior Arrests" , "Married","Violent","lsir Score","Age",
         "Custody Level","Prior Incarcerations",
         "Mental Health", "High School","Misconducts","Commit County")

sl.lib = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.earth","SL.xgboost","SL.glm")

#### with ranger classification for pi ####
# set up for train/test
set.seed(1234)
n = dim(df)[1]; nsplits = 2
s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

psis <- sds <- plugins <- pi.sds <- matrix(rep(NA,nsplits*4), ncol = nsplits)
assig.vecU <- assig.vecC<- assig.vecAr <- assig.vecAm <- assig.vecT <- infl.funU <- infl.funC <- infl.funAr <- infl.funAm <- rep(NA, n)
muFhat.mat <- matrix(rep(NA, n*length(avals)), ncol = length(avals))
pihat.old = read.csv('~jacquelinemauro/Dropbox/sorter/NoMediatorsSLpihat.csv')[,-1]
muhat.old = read.csv('~jacquelinemauro/Dropbox/sorter/NoMediatorsSLmuhat.csv')[,-1]
for(vfold in 1:nsplits){
  
  # make training/testing sets
  train = (s!=vfold); test = (s==vfold)
  train.df = df[train,]; test.df = df[test,]
  
  # output the constraint
  constr <- round(apply(pris.dummies[test,],2,sum) * 1.05)
  write.csv(constr, '~jacquelinemauro/Dropbox/sorter/capacityFudgeNm.csv')
  
  # get nuisance models and predictions
  mu.model = SuperLearner(Y = train.df$y, X = train.df[,-1], family = binomial(), SL.library = sl.lib)
  # newdat = test.df[,-1]
  # muhat <- matrix(rep(NA,length(avals)*dim(test.df)[1]), ncol = length(avals))
  # for(i in 1:length(avals)){
  #   newdat$A = avals[i]; newdat$A = factor(newdat$A, levels = avals)
  #   muhat[,i] = predict.SuperLearner(object = mu.model, newdata = newdat, onlySL = T)$pred
  # }
  pihat <- pihat.old[train,]
  muhat <- muhat.old[train,]
  # newdat = test.df[,-c(1,2)]; names(newdat) <- names(df)[-c(1,2)]
  # pi.model = ranger::ranger(as.factor(A)~., data = train.df[,-1], write.forest = T,probability = T)
  # pihat = predict(pi.model, newdat)$pre
  # 
  # # truncate extreme pi values
  # print(paste(sum(pihat < 1e-3 | pihat > 1-1e-3), 'extreme pi values'))
  # pihat[pihat < 1e-3] = 1e-3
  # pihat[pihat > 1-1e-3] = 1-1e-3
  # 
  # # output files
  # pihat.mat[test,] = pihat
  # muhat.mat[test,] = muhat
  # 
  # # get unconstrained plug in and if estimates
  # fU = apply(muhat,1,which.min)
  # pluginU = mean(apply(muhat,1,min))
  # fU.mat = sapply(c(1:length(avals)), function(a) as.numeric(fU == a))
  # pihatU = diag(pihat %*% t(fU.mat))
  # muhatU = diag(muhat %*% t(fU.mat))
  # ifU = (as.numeric(test.df$A == avals[fU])/pihatU)*(test.df$y - muhatU) + muhatU
  # psiU = mean(ifU)
  # sdU = sd(ifU)/sqrt(length(muhatU))
  
  # get constrained estimates
  # write.csv(constr, "~jacquelinemauro/Dropbox/sorter/capacityFudgeNm.csv")
  # write.csv(muhat, "~jacquelinemauro/Dropbox/sorter/SLmuhatUnconstrNewdatNmA.csv")
  # run_matlab_script("~jacquelinemauro/Dropbox/sorter/prison_assignment_sl_nm.m")
  
  # fC.mat = read.csv("~jacquelinemauro/Dropbox/sorter/prison_assignment_sl_nmA.csv", header = F)
  # fC = apply(fC.mat,1,which.max)
  # pihatC = diag(pihat %*% t(fC.mat))
  # muhatC = diag(muhat %*% t(fC.mat))
  # pluginC = mean(muhatC)
  # ifC = (as.numeric(test.df$A == avals[fC])/pihatC)*(test.df$y - muhatC) + muhatC
  # psiC = mean(ifC)
  # sdC = sd(ifC)/sqrt(length(muhatC))
  
  # get appromximate constrained regression-based estimates
  muFhat <- matrix(rep(NA,length(avals)*dim(train.df)[1]), ncol = length(avals))
  newdat = train.df[,-1]
  for(i in 1:length(avals)){
    newdat$A = avals[i]; newdat$A = factor(newdat$A, levels = avals)
    muFhat[,i] = predict.SuperLearner(object = mu.model, newdata = newdat, onlySL = T)$pred
  }
  muFhat.mat[train,] = muFhat
  
  constr.train <- round(apply(pris.dummies[train,],2,sum) * 1.05)
  write.csv(constr.train, "~jacquelinemauro/Dropbox/sorter/optConstr.csv")
  write.csv(muFhat, "~jacquelinemauro/Dropbox/sorter/optMuhat.csv")
  run_matlab_script("~jacquelinemauro/Dropbox/sorter/approxConstr.m")
  
  fAr.mat <- read.csv("~jacquelinemauro/Dropbox/sorter/optfhat.csv", header = F)
  fAr = apply(fAr.mat,1,which.max)
  fAr = factor(avals[fAr], levels = avals)
  
  class.df <- data.frame(A = fAr, muFhat)
  muF.model <- ranger::ranger(A~., data = class.df, write.forest = TRUE)
  fArhat <- predict(muF.model, data.frame(muhat), type='response')$pre
  fArhat.mat <- sapply(avals, function(a) as.numeric(fArhat == a))
  
  # pihatAr = diag(pihat %*% t(fArhat.mat))
  # muhatAr = diag(muhat %*% t(fArhat.mat))
  # pluginAr = mean(muhatAr)
  # ifAr = (as.numeric(test.df$A == avals[fArhat])/pihatAr)*(test.df$y - muhatAr) + muhatAr
  # psiAr = mean(ifAr)
  # sdAr = sd(ifAr)/sqrt(length(muhatAr))
  
  # get appromximate constrained matching-based estimates
  library(MatchIt)
  match.data <- data.frame(rbind(cbind(rep(0,dim(muFhat)[1]),muFhat),cbind(rep(1,dim(muhat)[1]),muhat)))
  names(match.data)[1]<-'Group'
  form = formula(paste('Group~', paste('X',c(2:(dim(muhat)[2]+1)), sep = "", collapse = "+"), sep = ""))
  match.it <- matchit(form, data = match.data, method="nearest", ratio=1)
  match.mat <- match.it$match.matrix
  fAm = fAr[as.numeric(match.mat)]
  
  fAmhat = factor(avals[fAm], levels = avals)
  fAmhat.mat <- sapply(avals, function(a) as.numeric(fAmhat == a))
  
  pihatAm = diag(pihat %*% t(fAmhat.mat))
  muhatAm = diag(muhat %*% t(fAmhat.mat))
  pluginAm = mean(muhatAm, na.rm = T)
  ifAm = (as.numeric(test.df$A == fAmhat)/pihatAm)*(test.df$y - muhatAm) + muhatAm
  psiAm = mean(ifAm, na.rm = T)
  sdAm = sd(ifAm, na.rm = T)/sqrt(length(muhatAm))
  
  # scalars
  psis[,vfold] = c(NA,NA,NA,psiAm)
  sds[,vfold] = c(NA,NA,NA,sdAm)
  
  plugins[,vfold] = c(NA,NA,NA,pluginAm)
  pi.sds[,vfold] = c(NA,NA,NA,sd(muhatAm)/sqrt(length(muhatAm)))
  
  # vectors
  # assig.vecU[test] = fU
  # assig.vecC[test] = fC
  # assig.vecAr[test] = fArhat
  assig.vecAm[test] = fAmhat
  # assig.vecT[test] = fAr
  # infl.funU[test] = ifU
  # infl.funC[test] = ifC
  # infl.funAr[test] = ifAr
  infl.funAm[test] = ifAm
  
}

if.results = cbind(apply(psis,1,mean),apply(sds,1,mean))
write.csv(if.results, 'NoMediatorsIF.csv')

pi.results = cbind(apply(plugins,1,mean),apply(pi.sds,1,mean))
write.csv(pi.results, 'NoMediatorsPI.csv')

assig.vecs = cbind(assig.vecU,assig.vecC,assig.vecAr,assig.vecAm, assig.vecT)
write.csv(assig.vecs, 'NoMediatorsassigvecs.csv')

write.csv(pihat.mat, 'NoMediatorspihat.csv')
write.csv(muhat.mat, 'NoMediatorsmuhat.csv')


##### redo with superlearner for pi ####
# set up for train/test
set.seed(1234)
n = dim(df)[1]; nsplits = 2
s <- sample(rep(1:nsplits,ceiling(n/nsplits))[1:n])

psis <- sds <- plugins <- pi.sds <- matrix(rep(NA,nsplits*4), ncol = nsplits)
assig.vecU <- assig.vecC<- assig.vecAr <- assig.vecAm <- assig.vecT <- infl.funU <- infl.funC <- infl.funAr <- infl.funAm <- rep(NA, n)
pihat.mat <- muhat.mat <- muFhat.mat <- matrix(rep(NA, n*length(avals)), ncol = length(avals))
for(vfold in 1:nsplits){
  
  # make training/testing sets
  train = (s!=vfold); test = (s==vfold)
  train.df = df[train,]; test.df = df[test,]
  
  # output the constraint
  constr <- round(apply(pris.dummies[test,],2,sum) * 1.05)
  write.csv(constr, '~jacquelinemauro/Dropbox/sorter/capacityFudgeNm.csv')
  
  # get nuisance models and predictions
  mu.model = SuperLearner(Y = train.df$y, X = train.df[,-1], family = binomial(), SL.library = sl.lib)
  newdat = test.df[,-1]
  muhat <- pihat <- matrix(rep(NA,length(avals)*dim(test.df)[1]), ncol = length(avals))
  for(i in 1:length(avals)){
    newdat$A = avals[i]; newdat$A = factor(newdat$A, levels = avals)
    muhat[,i] = predict.SuperLearner(object = mu.model, newdata = newdat, onlySL = T)$pred
  }
  newdat = test.df[,-c(1,2)]; names(newdat) <- names(df)[-c(1,2)]
  sl.lib.pi = c("SL.glmnet","SL.glm.interaction", "SL.mean","SL.ranger","SL.rpart","SL.lda","SL.xgboost","SL.polymars")
  for(i in 1:length(avals)){
    pi.model = SuperLearner(Y = as.numeric(train.df$A == avals[i]), X = train.df[,-c(1,2)], family = binomial(), SL.library = sl.lib.pi)
    pihat[,i] = predict.SuperLearner(object = pi.model, newdata = newdat, onlySL = T)$pred
  }
  
  # truncate extreme pi values
  print(paste(sum(pihat < 1e-3 | pihat > 1-1e-3), 'extreme pi values'))
  pihat[pihat < 1e-3] = 1e-3
  pihat[pihat > 1-1e-3] = 1-1e-3
  
  # output files
  pihat.mat[test,] = pihat
  muhat.mat[test,] = muhat
  
  # get unconstrained plug in and if estimates
  fU = apply(muhat,1,which.min)
  pluginU = mean(apply(muhat,1,min))
  fU.mat = sapply(c(1:length(avals)), function(a) as.numeric(fU == a))
  pihatU = diag(pihat %*% t(fU.mat))
  muhatU = diag(muhat %*% t(fU.mat))
  ifU = (as.numeric(test.df$A == avals[fU])/pihatU)*(test.df$y - muhatU) + muhatU
  psiU = mean(ifU)
  sdU = sd(ifU)/sqrt(length(muhatU))
  
  # get constrained estimates
  write.csv(constr, "~jacquelinemauro/Dropbox/sorter/capacityFudgeNm.csv")
  write.csv(muhat, "~jacquelinemauro/Dropbox/sorter/SLmuhatUnconstrNewdatNmA.csv")
  run_matlab_script("~jacquelinemauro/Dropbox/sorter/prison_assignment_sl_nm.m")
  
  fC.mat = read.csv("~jacquelinemauro/Dropbox/sorter/prison_assignment_sl_nmA.csv", header = F)
  fC = apply(fC.mat,1,which.max)
  pihatC = diag(pihat %*% t(fC.mat))
  muhatC = diag(muhat %*% t(fC.mat))
  pluginC = mean(muhatC)
  ifC = (as.numeric(test.df$A == avals[fC])/pihatC)*(test.df$y - muhatC) + muhatC
  psiC = mean(ifC)
  sdC = sd(ifC)/sqrt(length(muhatC))
  
  # get appromximate constrained regression-based estimates
  muFhat <- matrix(rep(NA,length(avals)*dim(train.df)[1]), ncol = length(avals))
  newdat = train.df[,-1]
  for(i in 1:length(avals)){
    newdat$A = avals[i]; newdat$A = factor(newdat$A, levels = avals)
    muFhat[,i] = predict.SuperLearner(object = mu.model, newdata = newdat, onlySL = T)$pred
  }
  muFhat.mat[train,] = muFhat
  
  constr.train <- round(apply(pris.dummies[train,],2,sum) * 1.05)
  write.csv(constr.train, "~jacquelinemauro/Dropbox/sorter/optConstr.csv")
  write.csv(muFhat, "~jacquelinemauro/Dropbox/sorter/optMuhat.csv")
  run_matlab_script("~jacquelinemauro/Dropbox/sorter/approxConstr.m")
  
  fAr.mat <- read.csv("~jacquelinemauro/Dropbox/sorter/optfhat.csv", header = F)
  fAr = apply(fAr.mat,1,which.max)
  fAr = factor(avals[fAr], levels = avals)
  
  class.df <- data.frame(A = fAr, muFhat)
  muF.model <- ranger::ranger(A~., data = class.df, write.forest = TRUE)
  fArhat <- predict(muF.model, data.frame(muhat), type='response')$pre
  fArhat.mat <- sapply(avals, function(a) as.numeric(fArhat == a))
  
  pihatAr = diag(pihat %*% t(fArhat.mat))
  muhatAr = diag(muhat %*% t(fArhat.mat))
  pluginAr = mean(muhatAr)
  ifAr = (as.numeric(test.df$A == avals[fArhat])/pihatAr)*(test.df$y - muhatAr) + muhatAr
  psiAr = mean(ifAr)
  sdAr = sd(ifAr)/sqrt(length(muhatAr))
  
  # get appromximate constrained matching-based estimates
  library(MatchIt)
  match.data <- data.frame(rbind(cbind(rep(0,dim(muFhat)[1]),muFhat),cbind(rep(1,dim(muhat)[1]),muhat)))
  names(match.data)[1]<-'Group'
  form = formula(paste('Group~', paste('X',c(2:(dim(muhat)[2]+1)), sep = "", collapse = "+"), sep = ""))
  match.it <- matchit(form, data = match.data, method="optimal", ratio=1)
  match.mat <- match.it$match.matrix
  fAm = fAr[as.numeric(match.mat)]
  
  fAmhat = factor(avals[fAm], levels = avals)
  fAmhat.mat <- sapply(avals, function(a) as.numeric(fAmhat == a))
  
  pihatAm = diag(pihat %*% t(fAmhat.mat))
  muhatAm = diag(muhat %*% t(fAmhat.mat))
  pluginAm = mean(muhatAm, na.rm = T)
  ifAm = (as.numeric(test.df$A == fAmhat)/pihatAm)*(test.df$y - muhatAm) + muhatAm
  psiAm = mean(ifAm, na.rm = T)
  sdAm = sd(ifAm, na.rm = T)/sqrt(length(muhatAm))
  
  # scalars
  psis[,vfold] = c(psiU,psiC,psiAr,psiAm)
  sds[,vfold] = c(sdU,sdC,sdAr,sdAm)
  
  plugins[,vfold] = c(pluginU,pluginC,pluginAr,pluginAm)
  pi.sds[,vfold] = c(sd(muhatU)/sqrt(length(muhatU)),sd(muhatC)/sqrt(length(muhatC)),
                     sd(muhatAr)/sqrt(length(muhatAr)),sd(muhatAm)/sqrt(length(muhatAm)))
  
  # vectors
  assig.vecU[test] = fU
  assig.vecC[test] = fC
  assig.vecAr[test] = fArhat
  assig.vecAm[test] = fAmhat
  assig.vecT[test] = fAr
  infl.funU[test] = ifU
  infl.funC[test] = ifC
  infl.funAr[test] = ifAr
  infl.funAm[test] = ifAm
  
}

if.results = cbind(apply(psis,1,mean),apply(sds,1,mean))
write.csv(if.results, '~jacquelinemauro/Dropbox/sorter/NoMediatorsSLifresults.csv')

pi.results = cbind(apply(plugins,1,mean),apply(pi.sds,1,mean))
write.csv(pi.results, '~jacquelinemauro/Dropbox/sorter/NoMediatorsSLpiresults.csv')

assig.vecs = cbind(assig.vecU,assig.vecC,assig.vecAr,assig.vecAm, assig.vecT)
write.csv(assig.vecs, '~jacquelinemauro/Dropbox/sorter/NoMediatorsSLassigvecs.csv')

write.csv(pihat.mat, '~jacquelinemauro/Dropbox/sorter/NoMediatorsSLpihat.csv')
write.csv(muhat.mat, '~jacquelinemauro/Dropbox/sorter/NoMediatorsSLmuhat.csv')

