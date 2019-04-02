# helper functions

plot.cace.double <- function(cace.object,...){
  library(ggplot2)
  df = data.frame(est = cace.object$Estimate, llp = cace.object$CI.lower.point,
                  ulp = cace.object$CI.upper.point, llu = cace.object$CI.lower.unif,
                  ulu = cace.object$CI.upper.unif, delta1 = cace.object$delta1,
                  delta2 = cace.object$delta2,
                  Name = paste("(",delta2,", ", delta1,")", sep = ""))
  p <- ggplot(df) + geom_line(aes(x = delta1, y = est, colour = "Estimate"))+
    geom_point(aes(x = delta1, y = est, colour = "Estimate"))+
    geom_ribbon(aes(ymin = llp, ymax = ulp, x = delta1, fill = "Pointwise"), alpha = .3) +
    geom_ribbon(aes(ymin = llu, ymax = ulu, x = delta1, fill = "Uniform"), alpha = .3) +
    geom_text(aes(x=delta1,y=est,label=Name),hjust=0, vjust=0, size = 3)+
    theme_bw()+
    geom_hline(yintercept = 0)+
    xlim(round(min(delta1*.9)), round(max(delta1*1.1)))+
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

simFunc <- function(N=5000,delta = 1, psi = 2, zmax = Inf, zmin = -Inf){
  y0 = rnorm(N,0,1); meanx = c(0,0,0,0)
  alpha = matrix(c(1,1,-1,-1),nrow = 4)
  x = matrix(unlist(lapply(meanx, function(x) rnorm(N,x,1))), nrow = N, byrow =T)
  z = rnorm(N, t(alpha)%*%t(x), sd = 2)
  pos.cond = (z + delta <= zmax); min.cond = (z - delta >= zmin)

  a = as.numeric(z >= y0)
  aplus = as.numeric(z + pos.cond*delta >= y0); amin = as.numeric(z - min.cond*delta >= y0)
  y = y0 + a*psi
  yplus = y0 + aplus*psi; ymin = y0 + amin*psi

  true.eff = psi
  true.ymean = psi*pnorm(z)
  true.amean = pnorm(z)
  true.ymean.plus = psi*pnorm(z+pos.cond*delta)
  true.amean.plus = pnorm(z+pos.cond*delta)
  true.ymean.min = psi*pnorm(z-min.cond*delta)
  true.amean.min = pnorm(z-min.cond*delta)
  true.z <- dnorm(z, mean = t(alpha)%*%t(x), sd = 2)
  true.z.min <- dnorm(z-min.cond*delta, mean = t(alpha)%*%t(x), sd = 2)
  true.z.plus <- dnorm(z+pos.cond*delta, mean = t(alpha)%*%t(x), sd = 2)

  return(data.frame(y,a,z,x,yplus,ymin,aplus,amin,
                    true.ymean,true.amean,true.ymean.plus,true.ymean.min,
                    true.amean.plus,true.amean.min,
                    true.z, true.z.min, true.z.plus))
}
