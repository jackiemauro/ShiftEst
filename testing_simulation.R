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

  #true.z[true.z < 1e-3] <- 1e-3; true.z[true.z.plus < 1e-3] <- 1e-3; true.z[true.z.min < 1e-3] <- 1e-3

  return(data.frame(y,a,z,x,yplus,ymin,aplus,amin,
                    true.ymean,true.amean,true.ymean.plus,true.ymean.min,
                    true.amean.plus,true.amean.min,
                    true.z, true.z.min, true.z.plus))
}

N = 5000 #eventually change to a vector
nsim = 100 #eventually change to 500
deltas = c(1,2,3) # eventually change to seq(from = 1, to = 5, length.out = 15)
error.mean = 1
K = c(1.99,3.99) # eventually change to c(1.99,2.99,3.99,5.99)
B = N^(1/K)

f.num <- function(df){
  mu0 = df$true.ymean #+ rnorm(N, 1, 1)/b
  mu.plus = df$true.ymean.plus + rnorm(N, 1, 1)/b
  mu.min = df$true.ymean.min + rnorm(N, 2, 1)/b

  ratio.plus = df$true.z.min/df$true.z + rnorm(N, 1, 1)/b
  ratio.plus[ratio.plus > 100] <- 100
  ratio.min = df$true.z.plus/df$true.z + rnorm(N, 2, 1)/b
  ratio.min[ratio.min > 100] <- 100

  return(ratio.plus*(df$y - mu0) + mu.plus - ratio.min*(df$y - mu0) - mu.min)
}
f.den <- function(df){
  la0 = df$true.amean #+ rnorm(N, 1, 1)/b
  la.plus = df$true.amean.plus + rnorm(N, 1, 1)/b
  la.min = df$true.amean.min + rnorm(N, 2, 1)/b

  ratio.plus = df$true.z.min/df$true.z + rnorm(N, 1, 1)/b
  ratio.plus[ratio.plus > 100] <- 100
  ratio.min = df$true.z.plus/df$true.z + rnorm(N, 2, 1)/b
  ratio.min[ratio.min > 100] <- 100

  return((ratio.plus - ratio.min)*(df$a - la0) + la.plus - la.min)
}
pi.num <- function(df){
  mu.plus = df$true.ymean.plus + rnorm(N, 1, 1)/b
  mu.min = df$true.ymean.min + rnorm(N, 2, 1)/b
  return(  mu.plus - mu.min )
}
pi.den <- function(df){
  la.plus = df$true.amean.plus + rnorm(N, 1, 1)/b
  la.min = df$true.amean.min + rnorm(N, 2, 1)/b
  return(  la.plus - la.min )
}

if.out <- pi.out <- data.frame(means = c(), sd = c(), lower = c(), upper = c(), rate = c())
for(rate in 1:length(K)){
  b = B[rate]
  plug.in <- if.est <- matrix(rep(NA, length(deltas)*nsim), ncol = length(deltas))
  for(i in 1:length(deltas)){
    data.list <- lapply( 1:nsim, function(x) simFunc(N=N, delta = deltas[i]) )
    plug.in[,i] <- unlist(lapply(data.list, function(x) mean(pi.num(x))/mean(pi.den(x))))
    if.est[,i] <- unlist(lapply(data.list, function(x) mean(f.num(x))/mean(f.den(x))))
  }

  pi.means = apply(plug.in, 2, mean); pi.sd = apply(plug.in,2,sd)
  pi.lower = pi.means - 1.96*pi.sd; pi.upper = pi.means + 1.96*pi.sd
  pi.out = rbind(pi.out, cbind(pi.means,pi.sd,pi.lower,pi.upper,rep(K[rate],length(deltas))))

  if.means = apply(if.est, 2, mean);if.sd = apply(if.est,2,sd)
  if.lower = if.means - 1.96*if.sd; if.upper = if.means + 1.96*if.sd
  if.out = rbind(if.out, cbind(if.means,if.sd,if.lower,if.upper,rep(K[rate],length(deltas))))
}


########### single shift works ##########

# make sure to clear each time
# going to add steps towards double shift and see if/when it breaks
# step 1: add no noise y.mean.min and no noise a.mean.min (as well as no noise ratio) -- works (sorta)
# step 2: add 0 mean noise to y.mean.min and a.mean.min -- also sorta works
# step 3: add N(1,1) noise to y.mean.min -- pulling away from truth at low values of delta
# step 4: add N(1,1) noise to a.mean.min -- fails badly
# step 5: returning to no noise y.mean.min and a.mean.min to make sure it still works -- still (kinda)
# step 6: add ratio noise N(2,1) -- still very low variance and pulling away at low delta
# step 7: add N(2,1) noise to mu0 & la0 terms in xi(-d) term -- breaks. if and pi almost identical, low variance and far at low deltas
# step 8: add N(1,1) noise to y.mean.min and a.mean.min -- no improvement
# TYPO in ratio for numerator, rerunning under same conditions as step8 -- variance better but if and pi still totally identical
# step 9: changing ratio error to N(1,1) -- slightly better but not good
# step 10: changing e.mean to 3 -- slightly better but not good (pi detaching more at r>3 but if still biased at low delta)
# step 11: N(-1,1) on muM and laM -- bad bad bad
# step 12: going back to zero noise on mu0 in xi(-d) term --works! (up to r=6 when it goes crazy)
# step 13: adding N(1,1) error on mu0 and la0 in xi(-d) term -- does badly at 4 and 3 but well at 6???
# step 14: error term has to be the same on mu0 and la0, changing to N(1,1) -- no good
# step 15: changing to mean 0 -- still bad
# step 16: no noise on any mu0/la0 term

# restarting
# no noise on anything -- works
# N(0,1) on mu0/la0; N(1,1) on muP/muM -- works
# add N(1,1) on laP/laM -- works (plug in does just as well on all these iterations)
# add N(1,1) on ratio terms -- works
# add N(2,1) on laP -- biased downwards
# setting laP noise back to N(1,1) and putting N(2,1) on muP -- biased up
# setting both laP and muP to N(2,1) --still biased down, IF and PI almost identical
# replacing mu0/la0 with truth in second term -- still biased down
# changing ratio errors and mu0/la0 to N(0,1), muM/laM to N(0,1), muP/laP to N(1,1)

# need ratP*mu0 to be close to muP. Know E(ratP*mu0) = muP
# ==> (ratP + rnorm(n,0,1)/B)*(mu0 + rnorm(n,0,1)) vs muP + rnorm(n,1,1)/B
# ==> muP + ratP*rnorm(n,0,1)/B + mu0*rnorm(n,0,1)/B + (rnorm(n,0,1)/B * rnorm(n,0,1)/B) vs muP + rnorm(n,1,1)/B


rm(list = ls())
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
  mu0 = df$true.ymean +rnorm(N,df$z,1)/B #+(rnorm(N,1,1)/B)
  muP = df$true.ymean.plus +rnorm(N,df$z+delta,1)/B#+(rnorm(N,2,1)/B)
  muM = df$true.ymean.min +rnorm(N,df$z-delta,1)/B#+(rnorm(N,1,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,1,1)/B)#+(rnorm(N,0,1)/B)
  ratP = df$true.z.plus/df$true.z +(rnorm(N,1,1)/B)#+(rnorm(N,0,1)/B)
  #(ratM*(df$y - (df$true.ymean+(rnorm(N,2,1)/B))) + muP) - (ratP*(df$y - mu0) + muM)
  (ratM*(df$y - mu0) + muP) - (ratP*(df$y - mu0) + muM)
}
f.den <- function(df){
  la0 = df$true.amean +rnorm(N,df$z,1)/B#+(rnorm(N,1,1)/B)
  laP = df$true.amean.plus +rnorm(N,df$z+delta,1)/B#+(rnorm(N,1,1)/B)
  laM = df$true.amean.min +rnorm(N,df$z-delta,1)/B#+(rnorm(N,1,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,1,1)/B)#+(rnorm(N,0,1)/B)
  ratP = df$true.z.plus/df$true.z+(rnorm(N,1,1)/B)#+(rnorm(N,0,1)/B)
  (ratM*(df$a - la0) + laP) - (ratP*(df$a - la0) + laM)
}
pi.num <- function(df){df$true.ymean.plus +rnorm(N,df$z+delta,1)/B#+(rnorm(N,2,1)/B)
  - df$true.ymean.min-rnorm(N,df$z-delta,1)/B#-(rnorm(N,1,1)/B)
  }
pi.den <- function(df){df$true.amean.plus+rnorm(N,df$z+delta,1)/B#+(rnorm(N,1,1)/B)
  - df$true.amean.min-rnorm(N,df$z-delta,1)/B#-(rnorm(N,1,1)/B)
  }
get_if <- function(df){mean(f.num(df))/mean(f.den(df))}
get_pi <- function(df){mean(pi.num(df))/mean(pi.den(df))}
single.delta <- function(delta){
  dfs = lapply(1:n.sim, function(x) simFunc(N=N, delta = delta))
  if.ests = unlist(lapply(dfs, function(df) get_if(df)))
  pi.ests = unlist(lapply(dfs, function(df) get_pi(df)))
  return(data.frame(if.ests = if.ests, pi.ests = pi.ests))
}

n.sim = 500
deltas = seq(.5,4.5,length.out = 15)
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

nsim = 100
K = c(1.99,2.99,3.99,5.99)
psi <- true.eff <- 2
N = 5000 # size of dataset
bootstrap.n = 10 # bootstrap samples
delta = seq(.5, 5.5, length = 15)
zmax = Inf; zmin = -Inf
e.mean <- 2 # mean of error term

pb <- txtProgressBar(min = 0, max = nsim*length(K), style = 3)
bigIFest <- bigIFsd <- bigPIest <- bigPIsd <- bigMBL <- bigMBU <- bigPWL <- bigPWU <- bigCover <- bigPWCover <- list()
bigMBLPI <- bigMBUPI <- bigPWLPI <- bigPWUPI <- bigCoverPI <- bigPWCoverPI <- bigMBL
for(j in 1:length(K)){
  r = K[j]
  print(paste('rate = ',r))
  B = N^(1/r)
  PI1 = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  IF1 = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  MBlower = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  MBupper <- PWlower <- PWupper <- ifvals <- pwcover <- pwcoverPI <- MBlower
  MBupperPI <- PWlowerPI <- PWupperPI <- MBlowerPI <- MBlower
  cover <- coverPI <- rep(NA,nsim)

  for(i in 1:nsim){
    #print(paste('i = ',i))
    datlist <- lapply(delta, function(d) simFunc(N, delta = d, psi = psi, zmax = zmax, zmin = zmin))

    PI1[i,] = unlist(lapply(datlist, function(d) mean(d$true.ymean.plus+(rnorm(N,2,1)/B) - (d$true.ymean.min+(rnorm(N,1,1)/B)))/
                              mean(d$true.amean.plus+(rnorm(N,1,1)/B) - (d$true.amean.min+(rnorm(N,1,1)/B)))))
    IF1[i,] = unlist(lapply(datlist, function(d) mean(f.num(d))/mean(f.den(d))))

    ll2 <- ul2 <- ul1 <- ll1 <- ll2PI <- ul2PI <- ll1PI <- ul1PI <- rep(NA,length(delta))
    for(jj in 1:length(delta)){
      psihat = IF1[i,jj]
      #ifvals = (f.num(datlist[[jj]]) - psihat*f.den(datlist[[jj]])) / mean(datlist[[jj]]$true.amean.plus - datlist[[jj]]$a)
      ifvals = (f.num(datlist[[jj]]) - psihat*f.den(datlist[[jj]])) / mean(datlist[[jj]]$true.amean.plus - datlist[[jj]]$true.amean)
      mb <- mult.bootstrap(est.eff=psihat,sigma=sd(ifvals),ifvals=ifvals,alpha=.05,n=N,nbs=bootstrap.n)
      ll2[jj] <- mb$ll2; ul2[jj] <- mb$ul2

      # pointwise confidence interval
      ll1[jj] <- psihat - qnorm(.975)*sd(ifvals)/sqrt(N)
      ul1[jj] <- psihat + qnorm(.975)*sd(ifvals)/sqrt(N)

      # corrected confidence intervals for plug in
      pi.est <- PI1[i,jj]
      ll1PI[jj] <- pi.est - qnorm(.975)*sd(ifvals)/sqrt(N)
      ul1PI[jj] <- pi.est + qnorm(.975)*sd(ifvals)/sqrt(N)
      ll2PI[jj] <- pi.est - (psihat - ll2[jj])
      ul2PI[jj] <- pi.est + (ul2[jj] - psihat)
    }

    MBlower[i,] = ll2; MBupper[i,] = ul2
    PWlower[i,] = ll1; PWupper[i,] = ul1

    MBlowerPI[i,] = ll2PI; MBupperPI[i,] = ul2PI
    PWlowerPI[i,] = ll1PI; PWupperPI[i,] = ul1PI

    cover[i] = (max(ll2) < psi) & (min(ul2) > psi)
    coverPI[i] = (max(ll2PI) < psi) & (min(ul2PI) > psi)
    pwcover[i,] = (ll1<psi) & (ul1>psi)
    pwcoverPI[i,] = (ll1PI<psi) & (ul1PI>psi)
    setTxtProgressBar(pb, (j-1)*nsim + i)
  }
  #write.csv(PI1, paste('tempPI_',j,sep = ""))
  #write.csv(IF1, paste('tempPI_',j,sep = ""))
  #write.csv(cbind(MBlower,MBupper), paste('tempMB_',j,sep = ""))
  #write.csv(cbind(MBlowerPI,MBupperPI), paste('tempMBPI_',j,sep = ""))

  bigPIest[[j]] = apply(PI1,2,mean); bigPIsd[[j]] = apply(PI1,2,sd)
  bigIFest[[j]] =  apply(IF1,2,mean); bigIFsd[[j]] = apply(IF1,2,sd)
  print(apply(IF1,2,mean))

  bigMBL[[j]] = MBlower; bigMBU[[j]] = MBupper
  bigPWL[[j]] = PWlower; bigPWU[[j]] = PWupper

  bigMBLPI[[j]] = MBlowerPI; bigMBUPI[[j]] = MBupperPI
  bigPWLPI[[j]] = PWlowerPI; bigPWUPI[[j]] = PWupperPI

  bigCover[[j]] = mean(cover)
  bigPWCover[[j]] = apply(pwcover,2,mean)

  bigCoverPI[[j]] = mean(coverPI)
  bigPWCoverPI[[j]] = apply(pwcoverPI,2,mean)
}
close(pb)



df = data.frame(upper.mb = rep(unlist(lapply(bigMBU, function(x) apply(x,2,mean))),2),
                lower.mb = rep(unlist(lapply(bigMBL, function(x) apply(x,2,mean))),2),
                upper.pw = rep(unlist(lapply(bigPWU, function(x) apply(x,2,mean))),2),
                lower.pw = rep(unlist(lapply(bigPWL, function(x) apply(x,2,mean))),2),
                est = c(unlist(bigPIest), unlist(bigIFest)),
                sd = c(unlist(bigPIsd), unlist(bigIFsd)),
                type = rep(c("PI","IF"), each = length(unlist(bigPIest))),
                delta = rep(delta, 2*length(K)),
                rate = rep(round(K,1), each = length(delta), times = 2)
)
df$lower.emp <- df$est - 1.96*df$sd
df$upper.emp <- df$est + 1.96*df$sd

emp.cover <- (df$lower.emp<2) & (df$upper.emp>2)

# coverage
cover <- data.frame(mb = rep(unlist(bigCover),each = length(delta)),
                    pw = unlist(bigPWCover),
                    emp.if = emp.cover[21:40],
                    emp.pi = emp.cover[1:20],
                    delta = rep(delta, length(K)),
                    rate = rep(round(K), each = length(delta)))
ggplot(df) +
  geom_hline(yintercept = psi, colour = 'red')+
  geom_point(aes(x = delta, y = est, shape = type))+
  geom_line(aes(x = delta, y = est, colour = type))+
  geom_ribbon(aes(ymin = lower.pw, ymax = upper.pw, x = delta, fill = "Pointwise"), alpha = .3) +
  geom_ribbon(aes(ymin = lower.mb, ymax = upper.mb, x = delta, fill = "Uniform"), alpha = .3)+
  facet_wrap(~rate)+
  theme_bw() +
  scale_color_manual(values = c('black', 'grey'))+
  theme(legend.position="bottom")+
  coord_cartesian(ylim = c(0, 3)) +
  labs(y = paste('Estimates (',nsim,' simulations)', sep = ""), x = "Shift Amount",
       title = "Estimates by estimator type, error rate and shift amount")

ggplot(df) +
  geom_hline(yintercept = psi, colour = 'red')+
  geom_point(aes(x = delta, y = est, shape = type))+
  geom_line(aes(x = delta, y = est, colour = type))+
  geom_ribbon(aes(ymin = lower.emp, ymax = upper.emp, x = delta, fill = type), alpha = .3) +
  facet_wrap(~rate)+
  theme_bw() +
  scale_color_manual(values = c('black', 'grey'))+
  scale_fill_manual(values = c('black', 'grey'))+
  theme(legend.position="bottom")+
  coord_cartesian(ylim = c(0, 3.5)) +
  labs(y = paste('Estimates (',nsim,' simulations)', sep = ""), x = "Shift Amount",
       title = "Estimates by estimator type, error rate and shift amount")
