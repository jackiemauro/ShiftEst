#' Single Shift Estimator Across Delta Range--obsolete, set delta2 = 0 in double shift
#'
#' @description Estimates the effects of shifting up 
#' across a range of values. Includes the multiplier bootstrap
#' for calculations of uniform confidence bands
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param x a dataframe of covariates
#' @param delta a vector of shift levels
#' @param Y.est an algorithm for estimating the means of your outcome.
#' Should be 'glm', 'superlearner' or 'ranger'. If you choose superlearner,
#' you can also specify the libraries you would like to use. Default libraries
#' are c("SL.glm","SL.randomForest","SL.polymars","SL.mean").
#' @param A.est an algorithm for estimating the means of your treatment.
#' Should be 'glm', 'superlearner' or 'ranger'. If you choose superlearner,
#' you can also specify the libraries you would like to use. Default libraries
#' are c("SL.glm","SL.randomForest","SL.polymars","SL.mean").
#' @param Z.est an algorithm for estimating the means of your treatment.
#' Should be 'glm' in which case a glm is used to estimte the mean and variance
#' and a kernel is used to estimate the density, or 'flexcode'. If you choose
#' 'flexcode', you can specify the regression function. The default is regressionFunction.NW
#' @param nfolds number of folds. Defaults to 2
#' @param zmax the upper bound on Z, default is Inf
#' @param zmin the lower bound on Z, default is -Inf
#' @param alpha the alpha level for your multiplier bootstrap confidence bands (default 0.05)
#' @param nbs the number of boostrap samples to take. Default is 10,000.
#' @param pos.cutoff the level at which you want to truncate pi(Z+delta)/pi(Z) for positivity
#' restrictions. Default is 50.
#'
#' @return a list including an estimate of the effect for each delta value and its standard deviation,
#' as well as upper and lower confidence intervals for the pointwise and uniform case.

# current issues:

single.shift.range <- function(y,a,z,x,delta,Y.est,A.est,Z.est,
                               nfolds = 2, zmax = Inf, zmin = -Inf,
                               alpha = 0.05, nbs = 10000,
                               pos.cutoff = 500, ...){
  
  #move these to their own function
  N = length(y); j = length(delta)
  if(any(length(a)!=N,length(z)!=N,dim(x)[1]!=N)){stop('y,a,z and x must be same length')}
  if(zmax <= zmin){stop('zmax must be bigger than zmin')}
  keep <- complete.cases(cbind(y,a,z,x))
  y <- y[keep]; a <- a[keep]; z <- z[keep]; x <- as.matrix(x[keep,])
  #need to drop out infinite values of x as well
  
  n = length(y)
  s = sample(rep(1:nfolds,ceiling(n/nfolds))[1:n])
  
  dat = as.data.frame(cbind(z, x))
  psihat <- sd <- compl <- matrix(rep(NA,nfolds*j), ncol = j, nrow = nfolds)
  ifvals <- matrix(rep(NA,n*j), ncol = j, nrow = n)
  
  for(i in 1:nfolds){
    train = (s!=i); test = (s==i)
    if(nfolds == 1){train <- test}
    
    ymean = y.mean.est(y[train],dat[train,],Y.est)
    amean = a.mean.est(a[train],dat[train,],A.est)
    zmean = z.condldens.est(z[train],x[train,],Z.est)
    
    mu.hat = my.pred(model = ymean, newdata = dat[test,], algo = Y.est)
    lambda.hat = my.pred(model = amean, newdata = dat[test,], algo = A.est)
    pi.hat = my.zpred(model = zmean, z = z[test], newdata = dat[test,-1], algo = Z.est)
    
    dat.plus = lapply(delta, function(d) as.data.frame(cbind(z+d,x)))
    dat.min = lapply(delta, function(d) as.data.frame(cbind(z-d,x)))
    
    for(ii in 1:length(dat.plus)){names(dat.plus[[ii]]) <- names(dat.min[[ii]]) <- names(dat)}
    
    mu.hat.plus <- lapply(dat.plus, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    lambda.hat.plus <- lapply(dat.plus, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))
    #pi.hat.plus <- lapply(delta, function(d) my.zpred(model = zmean, z = z[test]+d, newdata = dat[test,-1], algo = Z.est))
    
    #mu.hat.min <- lapply(dat.min, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    #lambda.hat.min <- lapply(dat.min, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))
    pi.hat.min <- lapply(delta, function(d) my.zpred(model = zmean, z = z[test]-d, newdata = dat[test,-1], algo = Z.est))
    
    for(jj in 1:j){
      d = delta[jj]
      mp = mu.hat.plus[[jj]];  #mm = mu.hat.min[[jj]]
      lp = lambda.hat.plus[[jj]];  #lm = lambda.hat.min[[jj]]
      #pp = pi.hat.plus[[jj]];  
      pm = pi.hat.min[[jj]]
      pos.cond = ((z[test]+d) < zmax)
      min.cond = ((z[test]-d) > zmin)
      
      if(sum( (abs(pm/pi.hat) > pos.cutoff), pi.hat==0, na.rm = T)>0){
        warning(print(paste(sum( (abs(pm/pi.hat) > pos.cutoff), pi.hat==0, na.rm = T),' extreme values found for pi ratio, setting to NA')))}
      
      pos.cond[(abs(pm/pi.hat) > pos.cutoff )] <- NA
      #min.cond[(abs(pp/pi.hat) > pos.cutoff )] <- NA
      
      xi.y.delta = ((y[test] - mu.hat)*(pm/pi.hat) - (y[test] - mp))*pos.cond
      xi.a.delta = ((a[test] - lambda.hat)*(pm/pi.hat) - (a[test] - lp))*pos.cond
      
      psihat[i,jj] = mean(xi.y.delta, na.rm = T)/mean(xi.a.delta, na.rm = T)
      ifvals[s==i,jj] = (xi.y.delta - psihat[i,jj]*(xi.a.delta))/mean(lp-lambda.hat, na.rm = T)
      compl[i,jj] = mean(xi.a.delta, na.rm = T)
      
    }
  }
  # average within delta values
  est.eff = apply(psihat, 2, mean, na.rm = T)
  sigma = sqrt(apply(ifvals,2, var, na.rm = T))
  compliers = apply(compl, 2, mean, na.rm = T)
  
  # pointwise confidence interval
  ll1 <- est.eff-qnorm(1-alpha/2)*sigma/sqrt(n); ul1 <- est.eff+qnorm(1-alpha/2)*sigma/sqrt(n)
  
  # multiplier bootstrap
  mb <- mult.bootstrap(est.eff=est.eff,sigma=sigma,ifvals=ifvals,alpha=alpha,n=n,nbs=nbs)
  calpha <- mb$calpha; ll2 <- mb$ll2; ul2 <- mb$ul2
  
  out = list(Estimate = est.eff, SD = sigma/sqrt(n), CI.lower.point = ll1, CI.upper.point = ul1,
             MBquantile = calpha, CI.lower.unif = ll2, CI.upper.unif = ul2, delta = delta, compliers = compliers)
  
  return(out)
}

mult.bootstrap <- function(est.eff,sigma,ifvals,alpha, n,nbs){
  eff.mat <- matrix(rep(est.eff,n), nrow = n, byrow = T)
  sig.mat <- matrix(rep(sigma,n), nrow = n, byrow = T)
  ifvals.std <- (ifvals - eff.mat)/sig.mat
  mult <- matrix(2*rbinom(n*nbs, 1, .5)-1, nrow = n, ncol = nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(apply(mult[,col]*ifvals.std,2,sum, na.rm = T)/sqrt(n)), na.rm = TRUE)
  })
  calpha <- quantile(maxvals, 1-alpha)
  ll2 <- est.eff - calpha*sigma/sqrt(n)
  ul2 <- est.eff + calpha*sigma/sqrt(n)
  
  return(list(calpha = calpha, ll2 = ll2, ul2 = ul2))
}

psi.est <- function(delta){
  dat.plus = as.data.frame(cbind(z+delta, x))
  dat.min = as.data.frame(cbind(z-delta, x))
  
  mu.hat.plus = my.pred(model = ymean, newdata = dat.plus[test,], algo = Y.est)
  lambda.hat.plus = my.pred(model = amean, newdata = dat.plus[test,], algo = A.est)
  pi.hat.plus = my.zpred(model = zmean, z = z[test], newdata = dat.plus[test,], algo = Z.est)
  
  mu.hat.min = my.pred(model = ymean, newdata = dat.min[test,], algo = Y.est)
  lambda.hat.min = my.pred(model = amean, newdata = dat.min[test,], algo = A.est)
  pi.hat.min = my.zpred(model = zmean, z = z[test], newdata = dat.min[test,], algo = Z.est)
  
  xi.y.delta = ((y[test] - mu.hat)*(pi.hat.min/pi.hat) - (y[test] - mu.hat.plus))*((z[test]+delta) < zmax)
  xi.y.Mdelta = ((y[test] - mu.hat)*(pi.hat.plus/pi.hat) - (y[test] - mu.hat.min))*((z[test]-delta) > zmin)
  
  xi.a.delta = ((a[test] - lambda.hat)*(pi.hat.min/pi.hat) - (a[test] - lambda.hat.plus))*((z[test]+delta) < zmax)
  xi.a.Mdelta = ((a[test] - lambda.hat)*(pi.hat.plus/pi.hat) - (a[test] - lambda.hat.min))*((z[test]-delta) > zmin)
  
  mean(xi.y.delta - xi.y.Mdelta)/mean(xi.a.delta - xi.a.Mdelta)
}

plot.cace <- function(cace.object,...){
  library(ggplot2)
  df = data.frame(est = cace.object$Estimate, llp = cace.object$CI.lower.point,
                  ulp = cace.object$CI.upper.point, llu = cace.object$CI.lower.unif,
                  ulu = cace.object$CI.upper.unif, delta = cace.object$delta)
  p <- ggplot(df) + geom_line(aes(x = delta, y = est, colour = "Estimate"))+
    geom_point(aes(x = delta, y = est, colour = "Estimate"))+
    geom_ribbon(aes(ymin = llp, ymax = ulp, x = delta, fill = "Pointwise"), alpha = .3) +
    geom_ribbon(aes(ymin = llu, ymax = ulu, x = delta, fill = "Uniform"), alpha = .3) +
    labs(x = expression(delta), y = expression(psi), title = "Confidence bands for Shift Estimator",
         fill = 'band type', colour = '')
  p
}


my.pred <- function(model, newdata, algo){
  if(algo == 'ranger'){predict(model, newdata)$pred}
  else if(algo == 'superlearner'){predict(model, newdata, onlySL = T)$pred}
  else{predict(model, newdata, type = 'response')}
}

my.zpred <- function(model, z, newdata, algo){
  if(algo == 'glm'){
    zhat = predict(model, newdata, type = 'response')
    z.var = mean( (z - zhat)^2  )
    N = length(zhat)
    
    gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}
    out = sapply(z, function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
  }
  else if(algo == 'ranger'){
    zhat = predict(model, newdata)$pred
    z.var = mean( (z - zhat)^2  )
    N = length(zhat)
    
    gK <- function(x){(1/sqrt(2*pi))*exp(-(x^2)/2)}
    out = sapply(z, function(y) (1/N)*sum(gK(sqrt( ((y - zhat))^2/z.var ) )))
  }
  else if(algo == 'flexcode'){
    pred = predict(model, newdata)
    out = get_probs(z, pred$z, pred$CDE)
  }
  # if(sum(out==0)>0){
  #   warning(print(paste(sum(out==0),' zero values found for pi, setting to NA')))
  #   out[out==0] <- NA}
  return(out)
}

get_probs<-function(z,zRange,mat){
  out = c(rep(NA,length(z)))
  for(i in 1:length(z)){
    r = rank(append(z[i],zRange),ties.method = 'last')[1]
    if(r>length(zRange)){out[i]=mat[i,(r-1)]}
    else if(r==1){out[i]=mat[i,1]}
    else{
      fl_val = mat[i,(r-1)]
      cl_val = mat[i,r]
      out[i] = mean(c(fl_val,cl_val))
    }
  }
  return(out)
}


