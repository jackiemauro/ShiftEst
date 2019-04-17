library(ggplot2)
source('R/doubleShiftRangeWrapper.R')
source('R/yMeanEst.R')
source('R/aMeanEst.R')
source('R/zCondlEst.R')

dat <- read.csv('~jacquelinemauro/MergedData.csv')[,-1]
dat$no.visits.last.loc <- (1 - dat$visitslastlocyn1)
dat$no.childvisits <- (1 - dat$childyn)
dat$maxtime <- apply(dat[,2:26],1,max)
dat <- dat[-which(dat$NCRecid.Event == 'reincarceration'),] #drop if they go back for parole violation
attach(dat)

startdummy = which(names(dat)=='alb'); enddummy = which(names(dat)=='pit')
options(na.action='na.pass')
county.f = factor(dat$CountyClass); county.dummies = model.matrix(~county.f)[,-c(1,2)]
minName.f = factor(dat$minTimeName); minName.dummies = model.matrix(~minName.f)[,-c(1,2)]

# using minimum distance prison as a dummy to indicate location
covs = cbind(white,loslastloc,minName.dummies,ageyrs,urban,priorarrests,married,violent,lsirscore,
             highschoolgrad,custody_level,numofpriorinc,mh,dat[,c(startdummy:(enddummy-1))],numoftotalmisconducts)

detach(dat)
old.names <- names(covs)
names(covs) <- sapply(c(1:dim(covs)[2]), function(k) paste('x',k,sep = ""))

#### summary stats ####
library(plyr)
library(xtable)

sumstats = ddply(dat, .(visitslastlocyn1), summarize,
                 recidivate = mean(NCRecid3,na.rm = T),
                 #minTime = mean(minTime,na.rm = T),
                 currentTime = mean(total_time,na.rm = T),
                 lengthofstay = mean(loslastloc,na.rm = T),
                 white = mean(white,na.rm = T),
                 age = mean(ageyrs,na.rm = T),
                 urban = mean(urban,na.rm = T),
                 married = mean(married,na.rm = T),
                 mentalhealth = mean(mh,na.rm = T),
                 highschool = mean(highschoolgrad,na.rm = T),
                 violent = mean(violent,na.rm = T),
                 score = mean(lsirscore,na.rm = T),
                 numberofvisits = mean(visitslastloc1,na.rm = T),
                 custodylevel  =mean(custody_level,na.rm = T),
                 priorarrests = mean(priorarrests,na.rm = T),
                 priorincarcerations = mean(numofpriorinc,na.rm = T),
                 misconducts = mean(numofmisconductslastloc,na.rm = T)
)

sds = ddply(dat, .(visitslastlocyn1), summarize,
            recidivate = sd(NCRecid3,na.rm = T),
            #minTime = sd(minTime,na.rm = T),
            currentTime = sd(total_time,na.rm = T),
            lengthofstay = sd(loslastloc,na.rm = T),
            white = sd(white,na.rm = T),
            age = sd(ageyrs,na.rm = T),
            urban = sd(urban,na.rm = T),
            married = sd(married,na.rm = T),
            mentalhealth = sd(mh,na.rm = T),
            highschool = sd(highschoolgrad,na.rm = T),
            violent = sd(violent,na.rm = T),
            score = sd(lsirscore,na.rm = T),
            numberofvisits = sd(visitslastloc1,na.rm = T),
            custodylevel  =sd(custody_level,na.rm = T),
            priorarrests = sd(priorarrests,na.rm = T),
            priorincarcerations = sd(numofpriorinc,na.rm = T),
            misconducts = sd(numofmisconductslastloc,na.rm = T)
)

out.tab = rbind(sumstats[1,],sds[1,],sumstats[2,],sds[2,])
library(xtable)
x = xtable(t(sumstats))
print(x, type = 'latex','summaryStatsNoparole.tex')
x = xtable(t(out.tab))
print(x, type = 'latex','summaryStatsNoparoleSds.tex')

# first, have delta1 = -delta2
delta1 <- c(20,40,60,90,120)
delta2 <- -delta1
output1 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                             Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 5, pos.cutoff = 100)
p <- plot.cace.double(output1)
write.csv(matrix(unlist(output1[-5]),ncol = length(delta1),byrow = T),file = 'dubshiftEven_5foldRR.csv')
ggsave(plot = p, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShiftEven_5foldRR.png',height = 4, width = 7)

# child visitation delta1 = -delta2
delta1 <- c(20,40,60,90,120)
delta2 <- -delta1
output1c <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.childvisits,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 5, pos.cutoff = 100)
p <- plot.cace.double(output1c)
write.csv(matrix(unlist(output1c[-5]),ncol = length(delta1),byrow = T),file = 'dubshiftEvenChild_5foldRR.csv')
ggsave(plot = p, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShiftEvenChild_5foldRR.png',height = 4, width = 7)



# not done:
# now do a moving 20 minute window
delta1 <- c(-100,-40,0,20,60,120)
delta2 <- c(-120,-60,-20,0,40,60)
output2 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 2,
                              pos.cutoff = 50)
p2 <- plot.cace.double(output2)
write.csv(matrix(unlist(output2[-5]),ncol = length(delta1),byrow = T),file = 'dubshiftWindow.csv')
ggsave(plot = p2, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShiftWindow.png',height = 4, width = 7)

# 30 minute window
delta1 <- c(-90,-60,-30,0,30,60,90,120)
delta2 <- c(-120,-90,-60,-30,0,30,60,90)
output3 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 2,
                              pos.cutoff = 50)
p3 <- plot.cace.double(output3)
write.csv(matrix(unlist(output3[-5]),ncol = length(delta1),byrow = T),file = 'dubshift30Window.csv')
ggsave(plot = p3, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShift30Window.png',height = 4, width = 7)

# 60 minute window
delta1 <- c(-120,-60,0,60,120,180)
delta2 <- c(-180,-120,-60,0,60,120)
output4 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 5,
                              pos.cutoff = 50)
p4 <- plot.cace.double(output4)
write.csv(matrix(unlist(output4[-5]),ncol = length(delta1),byrow = T),file = 'dubshift60Window.csv')
ggsave(plot = p4, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShift60Window.png',height = 4, width = 7)

############ run tsls #################
library(AER)
library(xtable)

df <- data.frame(subset(dat, select = c('NCRecid3', 'no.visits.last.loc', 'total_time')), covs)
df <- df[complete.cases(df),]
names(df) <- c('Recid', 'no.visits.last.loc', 'total_time', old.names)

# linear specification
tsls.form <- as.formula(paste("Recid~no.visits.last.loc+",paste(names(covs),collapse = "+"),"| total_time +",paste(names(covs),collapse = "+")))
tsls<- ivreg(tsls.form, data = df)
print(xtable(summary(tsls)$coefficients),type = "latex",file = "~/Dropbox/double robust causality/TSLS_output_RR.tex")

# more flexible specification -- just putting default splines on all non-binary variables
library(splines)
tsls.spline <- ivreg(Recid ~ no.visits.last.loc + white + bs(loslastloc) + minName.fCHS +
                       minName.fCOA + minName.fCRE + minName.fDAL + minName.fFRA +
                       minName.fFRS + minName.fFYT + minName.fGRA + minName.fGRE +
                       minName.fGRN + minName.fHOU + minName.fHUN + minName.fLAU +
                       minName.fMER + minName.fPIT + minName.fPNG + minName.fQUE +
                       minName.fRET + minName.fROC + minName.fSMR + minName.fWAM +
                       bs(ageyrs) + urban + bs(priorarrests) + married + violent + bs(lsirscore) +
                       highschoolgrad + custody_level + bs(numofpriorinc) + mh + alb +
                       cam + chs + coa + cre + dal + fra + frs + fyt + gra + gre +
                       grn + hou + hun + lau + mah + mer + png + que + ret + roc +
                       smi + smr + wam + bs(numoftotalmisconducts) | bs(total_time) + white +
                       bs(loslastloc) + minName.fCHS +
                       minName.fCOA + minName.fCRE + minName.fDAL + minName.fFRA +
                       minName.fFRS + minName.fFYT + minName.fGRA + minName.fGRE +
                       minName.fGRN + minName.fHOU + minName.fHUN + minName.fLAU +
                       minName.fMER + minName.fPIT + minName.fPNG + minName.fQUE +
                       minName.fRET + minName.fROC + minName.fSMR + minName.fWAM +
                       bs(ageyrs) + urban + bs(priorarrests) + married + violent + bs(lsirscore) +
                       highschoolgrad + custody_level + bs(numofpriorinc) + mh + alb +
                       cam + chs + coa + cre + dal + fra + frs + fyt + gra + gre +
                       grn + hou + hun + lau + mah + mer + png + que + ret + roc +
                       smi + smr + wam + bs(numoftotalmisconducts), data = df)

summary(tsls.spline)
print(xtable(summary(tsls.spline)$coefficients),type = "latex",file = "~/Dropbox/double robust causality/TSLS_splines_output_RR.tex")
