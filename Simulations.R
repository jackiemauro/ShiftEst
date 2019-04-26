rm(list = ls())

today = format(Sys.time(), "%Y%m%d")
set.seed(123)
library(ggplot2)
library(AER)

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
mult.bootstrap <- function(est.eff,sigma,ifvals,alpha, n,nbs){
  eff.mat <- matrix(rep(est.eff,n), nrow = n, byrow = T)
  sig.mat <- matrix(rep(sigma,n), nrow = n, byrow = T)
  ifvals.std <- (ifvals - eff.mat)/sig.mat
  mult <- matrix(2*rbinom(n*nbs, 1, .5)-1, nrow = n, ncol = nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(sum(mult[,col]*ifvals.std, na.rm = T)/sqrt(n)), na.rm = TRUE)
  })
  calpha <- quantile(maxvals, 1-alpha)
  ll2 <- est.eff - calpha*sigma/sqrt(n)
  ul2 <- est.eff + calpha*sigma/sqrt(n)

  return(list(calpha = calpha, ll2 = ll2, ul2 = ul2))
}
f.num <- function(df){
  mu0 = df$true.ymean +rnorm(N,df$z,1)/B
  muP = df$true.ymean.plus +rnorm(N,df$z + delta,1)/B
  muM = df$true.ymean.min +rnorm(N,df$z - delta,1)/B
  ratM = df$true.z.min/df$true.z +(rnorm(N,1,1)/B)
  ratP = df$true.z.plus/df$true.z +(rnorm(N,1,1)/B)
  (ratM*(df$y - mu0) + muP) - (ratP*(df$y - mu0) + muM)
}
f.den <- function(df){
  la0 = df$true.amean +rnorm(N,df$z,1)/B
  laP = df$true.amean.plus +rnorm(N,df$z + delta,1)/B
  laM = df$true.amean.min +rnorm(N,df$z - delta,1)/B
  ratM = df$true.z.min/df$true.z +(rnorm(N,1,1)/B)
  ratP = df$true.z.plus/df$true.z+(rnorm(N,1,1)/B)
  (ratM*(df$a - la0) + laP) - (ratP*(df$a - la0) + laM)
}
pi.num <- function(df){(df$true.ymean.plus+rnorm(N,df$z + delta,1)/B ) - (df$true.ymean.min +rnorm(N,df$z - delta,1)/B  )}
pi.den <- function(df){(df$true.amean.plus +rnorm(N,df$z + delta,1)/B) - (df$true.amean.min+rnorm(N,df$z - delta,1)/B)}
get_if <- function(df){mean(f.num(df))/mean(f.den(df))}
get_pi <- function(df){mean(pi.num(df))/mean(pi.den(df))}
single.delta <- function(delta){
  dfs = lapply(1:n.sim, function(x) simFunc(N=N, delta = delta))
  if.ests = unlist(lapply(dfs, function(df) get_if(df)))
  pi.ests = unlist(lapply(dfs, function(df) get_pi(df)))
  return(data.frame(if.ests = if.ests, pi.ests = pi.ests))
}

n.sim = 500
deltas = seq(.5,4,length.out = 15)
K = c(1.99,2.99,3.99,5.99)
Ns = c(100,1000,5000,10000)
psi <- true.eff <- 2
output <- list()
for(s.size in 1:length(Ns)){
  N = Ns[s.size]
  if.mean.by.rate <- pi.mean.by.rate <- if.sd.by.rate <- pi.sd.by.rate <- c()
  for(rate in 1:length(K)){
    print(round(K[rate]))
    B = N^(1/K[rate])
    if.ests <- pi.ests <- if.sds <- pi.sds <- c()
    for(d in 1:length(deltas)){
      delta = deltas[d]
      ests = single.delta(delta)
      if.ests[d] = mean(ests$if.ests); if.sds[d] = sd(ests$if.ests)
      pi.ests[d] = mean(ests$pi.ests); pi.sds[d] = sd(ests$pi.ests)
    }
    if.mean.by.rate = c(if.mean.by.rate, if.ests)
    pi.mean.by.rate = c(pi.mean.by.rate, pi.ests)
    if.sd.by.rate = c(if.sd.by.rate, if.sds)
    pi.sd.by.rate = c(pi.sd.by.rate, pi.sds)
  }
  if.out = data.frame(mean = if.mean.by.rate, sd = if.sd.by.rate, delta.val = rep(deltas, length(K)), rate.val = rep(round(K), each = length(deltas)))
  pi.out = data.frame(mean = pi.mean.by.rate, sd = pi.sd.by.rate, delta.val = rep(deltas, length(K)), rate.val = rep(round(K), each = length(deltas)))
  out = rbind(if.out, pi.out)
  out$type = c(rep("IF", dim(if.out)[1]),rep("PI", dim(pi.out)[1]))
  out$lower.emp = out$mean - 1.96*out$sd
  out$upper.emp = out$mean + 1.96*out$sd
  out$s.size = N

  output[[s.size]] = out
}

make.my.plot <- function(df){
  ggplot(df) +
    geom_hline(yintercept = psi, colour = 'red')+
    geom_point(aes(x = delta.val, y = mean, shape = type))+
    geom_line(aes(x = delta.val, y = mean, colour = type))+
    geom_ribbon(aes(ymin = lower.emp, ymax = upper.emp, x = delta.val, fill = type), alpha = .3) +
    facet_wrap(~rate.val)+
    theme_bw() +
    scale_color_manual(values = c('black', 'grey'))+
    scale_fill_manual(values = c('black', 'grey'))+
    theme(legend.position="bottom")+
    coord_cartesian(ylim = c(0, 4)) +
    labs(y = paste('Estimates (',n.sim,' simulations, sample size = ', df$s.size,')', sep = ""), x = "Shift Amount",
         title = "Estimates by estimator type, error rate and shift amount")
}
plots <- lapply(output, function(x) make.my.plot(x))
for(i in 1:length(plots)){ggsave(plot = plots[[i]], filename = paste("~/Dropbox/double robust causality/Figures/simulation_plot_",Ns[i],today,".png", sep = ""), height = 4, width = 7)}

############## april 25 edits ###########
# add closed form sd
# add tsls
# make sure it's over 500 sim's

rm(list = ls())

today = format(Sys.time(), "%Y%m%d")
set.seed(123)
library(ggplot2)
library(AER)


# need E(errorP(Z)) = E(error0(Z+delta)); E(errorM(Z)) = E(error0(Z-delta))
# this makes sense: error0(Z) = muhat(Z)-mu(Z) ==> error0(Z+d) = muhat(Z+d)-mu(Z+d) = errorP(Z)
error0 <- function(delta,df){rnorm(N, df$z, 1)/B}
errorP <- function(delta,df){rnorm(N,df$z+delta,1)/B}
errorM <- function(delta,df){rnorm(N,df$z-delta,1)/B}
rat.errorP <- function(delta,df){rnorm(N,1,1)/B}
rat.errorM <- function(delta,df){rnorm(N,1,1)/B}


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
mult.bootstrap <- function(est.eff,sigma,ifvals,alpha, n,nbs){
  eff.mat <- matrix(rep(est.eff,n), nrow = n, byrow = T)
  sig.mat <- matrix(rep(sigma,n), nrow = n, byrow = T)
  ifvals.std <- (ifvals - eff.mat)/sig.mat
  mult <- matrix(2*rbinom(n*nbs, 1, .5)-1, nrow = n, ncol = nbs)
  maxvals <- sapply(1:nbs, function(col){
    max(abs(sum(mult[,col]*ifvals.std, na.rm = T)/sqrt(n)), na.rm = TRUE)
  })
  calpha <- quantile(maxvals, 1-alpha)
  ll2 <- est.eff - calpha*sigma/sqrt(n)
  ul2 <- est.eff + calpha*sigma/sqrt(n)

  return(list(calpha = calpha, ll2 = ll2, ul2 = ul2))
}
f.num <- function(df){
  mu0 = df$true.ymean + error0(delta,df)
  muP = df$true.ymean.plus + errorP(delta,df)
  muM = df$true.ymean.min + errorM(delta,df)
  ratM = df$true.z.min/df$true.z + rat.errorM(delta,df)
  ratP = df$true.z.plus/df$true.z + rat.errorP(delta,df)
  (ratM*(df$y - mu0) + muP) - (ratP*(df$y - mu0) + muM)
}
f.den <- function(df){
  la0 = df$true.amean + error0(delta,df)
  laP = df$true.amean.plus + errorP(delta,df)
  laM = df$true.amean.min + errorM(delta,df)
  ratM = df$true.z.min/df$true.z + rat.errorM(delta,df)
  ratP = df$true.z.plus/df$true.z + rat.errorP(delta,df)
  (ratM*(df$a - la0) + laP) - (ratP*(df$a - la0) + laM)
}
pi.num <- function(df){(df$true.ymean.plus+errorP(delta,df)) - (df$true.ymean.min + errorM(delta,df))}
pi.den <- function(df){(df$true.amean.plus + errorP(delta,df) ) - (df$true.amean.min + errorM(delta,df))}
get_if <- function(df){mean(f.num(df))/mean(f.den(df))}
get_pi <- function(df){mean(pi.num(df))/mean(pi.den(df))}
get_tsls_est <- function(df){summary(ivreg(y~a+X1+X2+X3+X3 | z+X1+X2+X3+X3, data = df))$coef["a","Estimate"]}
get_tsls_sd <- function(df){summary(ivreg(y~a+X1+X2+X3+X3 | z+X1+X2+X3+X3, data = df))$coef["a","Std. Error"]}
psi_sd <- function(df){sd((f.num(df) - get_if(df)*f.den(df))/mean((df$true.amean.plus+errorM(delta,df))-(df$true.amean.min+errorP(delta,df))))/sqrt(dim(df)[1])}
single.delta <- function(delta){
  dfs = lapply(1:n.sim, function(x) simFunc(N=N, delta = delta))
  if.ests = unlist(lapply(dfs, function(df) get_if(df)))
  pi.ests = unlist(lapply(dfs, function(df) get_pi(df)))
  th.sds = unlist(lapply(dfs, function(df) psi_sd(df)))
  tsls.ests = unlist(lapply(dfs, function(df) get_tsls_est(df)))
  tsls.sds = unlist(lapply(dfs, function(df) get_tsls_sd(df)))
  return(data.frame(if.ests = if.ests, pi.ests = pi.ests, th.sds = th.sds, tsls.ests = tsls.ests, tsls.sds = tsls.sds))
}

n.sim = 500
deltas = seq(.5,4,length.out = 15)
K = c(1.99,2.99,3.99,5.99)
Ns = c(100,1000,5000,10000)

psi <- true.eff <- 2
output <- list()
for(s.size in 1:length(Ns)){
  N = Ns[s.size]
  if.mean.by.rate <- pi.mean.by.rate <- if.sd.by.rate <- pi.sd.by.rate <- th.sd.by.rate <- tsls.sd.by.rate <- tsls.mean.by.rate <- c()
  for(rate in 1:length(K)){
    print(paste("sample size = ",N,", rate = ",round(K[rate])), sep = "")
    B = N^(1/K[rate])
    if.ests <- pi.ests <- if.sds <- pi.sds <- th.sds <- tsls.ests <- tsls.sds <- c()
    for(d in 1:length(deltas)){
      delta = deltas[d]
      ests = single.delta(delta)
      th.sds[d] = mean(ests$th.sds)
      tsls.ests[d] = mean(ests$tsls.ests); tsls.sds[d] = mean(ests$tsls.ests)
      if.ests[d] = mean(ests$if.ests); if.sds[d] = sd(ests$if.ests)
      pi.ests[d] = mean(ests$pi.ests); pi.sds[d] = sd(ests$pi.ests)
    }
    if.mean.by.rate = c(if.mean.by.rate, if.ests)
    pi.mean.by.rate = c(pi.mean.by.rate, pi.ests)
    tsls.mean.by.rate = c(tsls.mean.by.rate, tsls.ests)
    if.sd.by.rate = c(if.sd.by.rate, if.sds)
    th.sd.by.rate = c(th.sd.by.rate, th.sds)
    pi.sd.by.rate = c(pi.sd.by.rate, pi.sds)
    tsls.sd.by.rate = c(tsls.sd.by.rate, tsls.sds)
  }
  if.out = data.frame(mean = if.mean.by.rate, sd = if.sd.by.rate, th.sd = th.sd.by.rate, delta.val = rep(deltas, length(K)), rate.val = rep(round(K), each = length(deltas)))
  pi.out = data.frame(mean = pi.mean.by.rate, sd = pi.sd.by.rate, th.sd = th.sd.by.rate, delta.val = rep(deltas, length(K)), rate.val = rep(round(K), each = length(deltas)))
  tsls.out = data.frame(mean = tsls.mean.by.rate, sd = tsls.sd.by.rate, th.sd = th.sd.by.rate, delta.val = rep(deltas, length(K)), rate.val = rep(round(K), each = length(deltas)))
  out = rbind(if.out, pi.out, tsls.out)
  out$type = c(rep("IF", dim(if.out)[1]),rep("PI", dim(pi.out)[1]),rep("TSLS", dim(tsls.out)[1]))
  out$lower.th = out$mean - 1.96*out$th.sd
  out$upper.th = out$mean + 1.96*out$th.sd
  out$lower.emp = out$mean - 1.96*out$sd
  out$upper.emp = out$mean + 1.96*out$sd
  out$s.size = N

  output[[s.size]] = out
}

make.my.plot <- function(df){
  ggplot(df) +
    geom_hline(yintercept = psi, colour = 'red')+
    geom_point(aes(x = delta.val, y = mean, shape = type))+
    geom_line(aes(x = delta.val, y = mean, colour = type))+
    geom_ribbon(aes(ymin = lower.emp, ymax = upper.emp, x = delta.val, fill = type), alpha = .3) +
    #geom_ribbon(aes(ymin = lower.th, ymax = upper.th, x = delta.val, col = type), alpha = .1) +
    facet_wrap(~rate.val)+
    theme_bw() +
    scale_color_manual(values = c('black', 'grey'))+
    scale_fill_manual(values = c('black', 'grey'))+
    theme(legend.position="bottom")+
    coord_cartesian(ylim = c(0, 4)) +
    labs(y = paste('Estimates (',n.sim,' simulations, sample size = ', df$s.size,')', sep = ""), x = "Shift Amount",
         title = "Estimates by estimator type, error rate and shift amount")
}
plots <- lapply(output, function(x) make.my.plot(x[x$type != "TSLS",]))
for(i in 1:length(plots)){ggsave(plot = plots[[i]], filename = paste("~/Dropbox/double robust causality/Figures/simulation_plot_",Ns[i],today,".png", sep = ""), height = 4, width = 7)}
lapply(output, function(k) write.csv(k,file=paste("~/Dropbox/double robust causality/df_",k$s.size[1],"_",today,".csv",sep="")))

