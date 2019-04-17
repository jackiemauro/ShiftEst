#' Double Shift Estimator Across Delta Range -- not right yet
#'
#' @description Estimates the effects of two shift amounts
#' across a range of values. Includes the multiplier bootstrap
#' for calculations of uniform confidence bands
#'
#' @param y your outcome variable (a vector)
#' @param a your treatment variable (a vector)
#' @param z your instrument (a vector)
#' @param x a dataframe of covariates
#' @param delta1 a vector of shift levels
#' @param delta2 a vector of shift levels, should be smaller than delta1.
#' If delta2 is not specified, it will automatically set to be -delta1. For single
#' shift, should be able to set this to 0.
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
# note that this is E(Y^1 - Y^0 | A^{Z+d1} > A^{Z+d2}) rather than Z+delta > Z-delta
# if delta2 is not specified, we take delta2 = -delta1 to recover Z+delta > Z-delta
#need to drop out infinite values of x as well

double.shift.range <- function(y,a,z,x,delta1,delta2=NULL,Y.est,A.est,Z.est,
                               nfolds = 2, zmax = Inf, zmin = -Inf,
                               alpha = 0.05, nbs = 10000,
                               pos.cutoff = 500, ...){

  # initial data checks
  N = length(y); j = length(delta1)
  if(any(length(a)!=N,length(z)!=N,dim(x)[1]!=N)){stop('y,a,z and x must be same length')}
  if(zmax <= zmin){stop('zmax must be bigger than zmin')}
  keep <- complete.cases(cbind(y,a,z,x))
  y <- y[keep]; a <- a[keep]; z <- z[keep]; x <- as.matrix(x[keep,])
  if(is.null(delta2)){delta2 <- (-delta1)}
  if(length(delta1)!=length(delta2)){stop('delta vectors must be same length')}
  if(any(is.infinite(x))){stop('infinite values detected in x')}

  # basic parameters
  n = length(y)
  s = sample(rep(1:nfolds,ceiling(n/nfolds))[1:n])
  dat = as.data.frame(cbind(z, x))

  # get parameter across folds
  psihat <- sd <- compl <- matrix(rep(NA,nfolds*j), ncol = j, nrow = nfolds)
  ifvals <- matrix(rep(NA,n*j), ncol = j, nrow = n)
  for(i in 1:nfolds){
    train = (s!=i); test = (s==i)
    if(nfolds == 1){train <- test}

    # train models
    ymean = y.mean.est(y[train],dat[train,],Y.est)
    amean = a.mean.est(a[train],dat[train,],A.est)
    zmean = z.condldens.est(z[train],x[train,],Z.est)

    # predict models on test data at observed z
    mu.hat = my.pred(model = ymean, newdata = dat[test,], algo = Y.est)
    lambda.hat = my.pred(model = amean, newdata = dat[test,], algo = A.est)
    pi.hat = my.zpred(model = zmean, z = z[test], newdata = dat[test,-1], algo = Z.est)

    # shifted datasets
    dat1 = lapply(delta1, function(d) as.data.frame(cbind(z+d,x)))
    dat2 = lapply(delta2, function(d) as.data.frame(cbind(z+d,x)))
    for(ii in 1:length(dat1)){names(dat1[[ii]]) <- names(dat2[[ii]]) <- names(dat)}

    # predict models on test data at delta1 shift
    mu1 <- lapply(dat1, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    lambda1 <- lapply(dat1, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))
    pi1 <- lapply(delta1, function(d) my.zpred(model = zmean, z = z[test]-d, newdata = dat[test,-1], algo = Z.est))

    # predict models on test data at delta2 shift
    mu2 <- lapply(dat2, function(df) my.pred(model = ymean, newdata = df[test,], algo = Y.est))
    lambda2 <- lapply(dat2, function(df) my.pred(model = amean, newdata = df[test,], algo = A.est))
    pi2 <- lapply(delta2, function(d) my.zpred(model = zmean, z = z[test]-d, newdata = dat[test,-1], algo = Z.est))

    for(jj in 1:j){
      d1 = delta1[jj]; d2 = delta2[jj]
      m1 = mu1[[jj]]; m2 = mu2[[jj]]
      l1 = lambda1[[jj]]; l2 = lambda2[[jj]]
      p1 = pi1[[jj]]; p2 = pi2[[jj]]

      pos.cond1 = (z[test]+d1<=zmax) & (z[test]+d1>=zmin)
      pos.cond2 = (z[test]+d2<=zmax) & (z[test]+d2>=zmin)

      if(sum( (abs(p2/pi.hat) > pos.cutoff), (abs(p1/pi.hat) > pos.cutoff), pi.hat==0, na.rm = T)>0){
        warning(print(paste(sum( (abs(p1/pi.hat) > pos.cutoff), (abs(p1/pi.hat) > pos.cutoff), pi.hat==0, na.rm = T),' extreme values found for pi ratio, setting to NA')))}
      pos.cond1[(abs(p1/pi.hat) > pos.cutoff )] <- NA
      pos.cond2[(abs(p2/pi.hat) > pos.cutoff )] <- NA

      xi1.y.delta = ((y[test] - mu.hat)*(p1/pi.hat) - (y[test]-m1) )*pos.cond1
      xi1.a.delta = ((a[test] - lambda.hat)*(p1/pi.hat) -(a[test]-l1))*pos.cond1
      xi2.y.delta = ((y[test] - mu.hat)*(p2/pi.hat) - (y[test]-m2))*pos.cond2
      xi2.a.delta = ((a[test] - lambda.hat)*(p2/pi.hat) - (a[test]-l2))*pos.cond2

      xi.y.delta = xi1.y.delta - xi2.y.delta
      xi.a.delta = xi1.a.delta - xi2.a.delta

      psihat[i,jj] = mean(xi.y.delta, na.rm = T)/mean(xi.a.delta, na.rm = T)
      ifvals[s==i,jj] = (xi.y.delta - psihat[i,jj]*(xi.a.delta))/mean(l1-l2, na.rm = T)
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
             MBquantile = calpha, CI.lower.unif = ll2, CI.upper.unif = ul2,
             delta1 = delta1, delta2 = delta2, compliers = compliers
  )

  return(out)
}

