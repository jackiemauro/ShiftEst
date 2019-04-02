rm(list = ls())
library(ggplot2)

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
  mu0 = df$true.ymean+(rnorm(N,e.mean,1)/B) 
  muP = df$true.ymean.plus+(rnorm(N,e.mean,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean,1)/B)
  ratM*(df$y - mu0) + muP
}
f.den <- function(df){
  la0 = df$true.amean+(rnorm(N,e.mean,1)/B) 
  laP = df$true.amean.plus+(rnorm(N,e.mean,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean,1)/B)
  ratM*(df$a - la0) + laP
}

### run some simulations adding various noise ----
nsim = 100
K = c(1.99,2.99,3.99,5.99)
psi <- true.eff <- 2
N = 5000 # size of dataset
bootstrap.n = 1e4 # bootstrap samples
delta = seq(.5, 5.5, length = 15)
zmax = Inf; zmin = -Inf
e.mean <- 2 # mean of error term

pb <- txtProgressBar(min = 0, max = nsim*length(K), style = 3)
bigIFest <- bigIFsd <- bigPIest <- bigPIsd <- bigMBL <- bigMBU <- bigPWL <- bigPWU <- bigCover <- bigPWCover <- list()
bigMBLPI <- bigMBUPI <- bigPWLPI <- bigPWUPI <- bigCoverPI <- bigPWCoverPI <- bigMBL
for(j in 1:length(K)){
  r = K[j]
  #print(paste('rate = ',r))
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
    
    PI1[i,] = unlist(lapply(datlist, function(d) mean(d$true.ymean.plus+(rnorm(N,e.mean,1)/B) - d$y)/mean(d$true.amean.plus+(rnorm(N,e.mean,1)/B) - d$a)))
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
  write.csv(PI1, paste('tempPI_',j,sep = ""))
  write.csv(IF1, paste('tempPI_',j,sep = ""))
  write.csv(cbind(MBlower,MBupper), paste('tempMB_',j,sep = ""))
  write.csv(cbind(MBlowerPI,MBupperPI), paste('tempMBPI_',j,sep = ""))

  bigPIest[[j]] = apply(PI1,2,mean); bigPIsd[[j]] = apply(PI1,2,sd)
  bigIFest[[j]] =  apply(IF1,2,mean); bigIFsd[[j]] = apply(IF1,2,sd)
  
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
write.csv(cover, file = 'simulationCoveragemu2.csv')

#write.csv(df, file = 'simulationOutput.csv') # mu = 1
write.csv(df, file = 'simulationOutputmu2.csv') # cranked error mu to 2
g<-ggplot(df) + 
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
#ggsave('type_rate_shift_plots.png',g,height = 4, width = 7) #mu=1
ggsave('type_rate_shift_plotsmu2.png',g,height = 4, width = 7) #mu=2

g<-ggplot(df) + 
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
#ggsave('type_rate_shift_plots_empSD.png',g,height = 4, width = 7)
ggsave('type_rate_shift_plots_empSDmu2.png',g,height = 4, width = 7)

#### figures at dense deltas
df2<- data.frame(upper.mb = c(unlist(lapply(bigMBUPI, function(x) apply(x,2,mean))),unlist(lapply(bigMBU, function(x) apply(x,2,mean)))),
                   lower.mb = c(unlist(lapply(bigMBLPI, function(x) apply(x,2,mean))),unlist(lapply(bigMBL, function(x) apply(x,2,mean)))),
                   upper.pw = c(unlist(lapply(bigPWUPI, function(x) apply(x,2,mean))),unlist(lapply(bigPWU, function(x) apply(x,2,mean)))),
                   lower.pw = c(unlist(lapply(bigPWLPI, function(x) apply(x,2,mean))),unlist(lapply(bigPWL, function(x) apply(x,2,mean)))),
                   est = c(unlist(bigPIest), unlist(bigIFest)),
                   sd = c(unlist(bigPIsd), unlist(bigIFsd)),
                   type = rep(c("PI","IF"), each = length(unlist(bigPIest))),
                   delta = rep(delta, 2*length(K)),
                   rate = rep(round(K,1), each = length(delta), times = 2)
)
write.csv(df2, file = 'simulationOutputmu2dense.csv') # cranked error mu to 2

g<-ggplot(df2) + 
  geom_hline(yintercept = psi, colour = 'red')+
  geom_line(aes(x = delta, y = est, colour = type))+
  geom_point(aes(x = delta, y = est, shape = type))+
  geom_ribbon(aes(ymin = lower.mb, ymax = upper.mb, x = delta, fill = type), alpha = .3) +
  facet_wrap(~rate)+ 
  theme_bw() + 
  scale_color_manual(values = c('grey', 'black'))+
  scale_fill_manual(values = c('grey', 'black'))+
  theme(legend.position="bottom")+
  coord_cartesian(ylim = c(0, 3.5)) +
  labs(y = paste('Estimates (',nsim,' simulations)', sep = ""), x = "Shift Amount", 
       title = "Estimates by estimator type, error rate and shift amount")
ggsave('type_rate_shift_plots_empSDmu2dense.png',g,height = 4, width = 7)
#### dense delta vs bias plot ----
rm(list = ls())
library(plyr)
library(ggplot2)

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
f.num <- function(df,B){
  mu0 = df$true.ymean+(rnorm(N,e.mean,1)/B) 
  muP = df$true.ymean.plus+(rnorm(N,e.mean,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean,1)/B)
  ratM*(df$y - mu0) - (df$y - muP)
}
f.den <- function(df,B){
  la0 = df$true.amean+(rnorm(N,e.mean,1)/B) 
  laP = df$true.amean.plus+(rnorm(N,e.mean,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean,1)/B)
  ratM*(df$a - la0) - (df$a - laP)
}

nsim = 100
K = c(1.99,2.99,3.99,5.99)
psi <- true.eff <- 2
N = 5000 # size of dataset
zmax = Inf; zmin = -Inf
e.mean <- 1 # mean of error term
delta = seq(1,10, length = 75)

temp <- matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
means <- sds <- bias <- abs.bias <- matrix(rep(NA,length(K)*length(delta)), ncol = length(delta))
for(B in 1:length(K)){
  for(i in 1:nsim){
    datlist <- lapply(delta, function(d) simFunc(N, delta = d, psi = psi, zmax = zmax, zmin = zmin))  
    temp[i,] = unlist(lapply(datlist, function(d) mean(f.num(d,B=N^(1/K[B])))/mean(f.den(d,B = N^(1/K[B])))))
  }
  bias[B,] = apply(temp-psi,2,mean,na.rm = T)
  abs.bias[B,] = apply(abs(temp-psi),2,mean,na.rm = T)
  means[B,] = apply(temp,2,mean,na.rm = T)
  sds[B,] = apply(temp,2,sd,na.rm = T)
}

df = data.frame(means = as.numeric(t(means)),
                sds = as.numeric(t(sds)),
                bias = as.numeric(t(bias)),
                abs.bias = as.numeric(t(abs.bias)),
                rate = rep(as.factor(round(K,1)), each = length(delta)),
                delta = rep(delta, length(K))
                )


g <- ggplot(df) + 
  geom_line(aes(x = delta, y = means, color = rate)) + 
  theme_bw() +
  ylim(c(-10,10)) +
  geom_hline(yintercept = psi)

ggsave(filename = 'rates_v_delta.png', g, height = 4, width = 7)

g <- ggplot(df) + 
  geom_line(aes(x = delta, y = bias, color = rate)) + 
  theme_bw() +
  ylim(c(-10,10)) +
  geom_hline(yintercept = 0)

ggsave(filename = 'bias_v_delta.png', g, height = 4, width = 7)

g <- ggplot(df) + 
  geom_line(aes(x = delta, y = abs.bias, color = rate)) + 
  theme_bw() +
  ylim(c(-10,10)) +
  geom_hline(yintercept = 0)

ggsave(filename = 'absbias_v_delta.png', g, height = 4, width = 7)
### run some simulations of double shift adding assym noise ----
f.num2 <- function(df){
  mu0 = df$true.ymean+(rnorm(N,e.mean1,1)/B) 
  muP = df$true.ymean.plus+(rnorm(N,e.mean1,1)/B)
  muM = df$true.ymean.min+(rnorm(N,e.mean2,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean2,1)/B)
  ratP = df$true.z.plus/df$true.z +(rnorm(N,e.mean1,1)/B)
  ratM*(df$y - mu0) + muP - ratP*(df$y - mu0) - muM
}
f.den2 <- function(df){
  la0 = df$true.amean+(rnorm(N,e.mean1,1)/B) 
  laP = df$true.amean.plus+(rnorm(N,e.mean1,1)/B)
  laM = df$true.amean.min+(rnorm(N,e.mean2,1)/B)
  ratM = df$true.z.min/df$true.z +(rnorm(N,e.mean2,1)/B)
  ratP = df$true.z.plus/df$true.z +(rnorm(N,e.mean1,1)/B)
  ratM*(df$a - la0) + laP - ratP*(df$a - la0) - laM
}

nsim = 100
K = c(1.99,2.99,3.99,5.99)
psi <- true.eff <- 2
N = 5000 # size of dataset
bootstrap.n = 10 # bootstrap samples
delta = seq(1, 5, length = 15)
zmax = Inf; zmin = -Inf
# for double shift, can't add symmetric error or it'll cancel largely
e.mean1 = 2
e.mean2 = 1

pb <- txtProgressBar(min = 0, max = nsim*length(K), style = 3)
bigIFest <- bigIFsd <- bigPIest <- bigPIsd <- bigMBL <- bigMBU <- bigPWL <- bigPWU <- bigCover <- bigPWCover <- list()
bigMBLPI <- bigMBUPI <- bigPWLPI <- bigPWUPI <- bigCoverPI <- bigPWCoverPI <- bigMBL
for(j in 1:length(K)){
  r = K[j]
  B = N^(1/r)
  PI1 = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  IF1 = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  MBlower = matrix(rep(NA,nsim*length(delta)), ncol = length(delta))
  MBupper <- PWlower <- PWupper <- ifvals <- pwcover <- pwcoverPI <- MBlower
  MBupperPI <- PWlowerPI <- PWupperPI <- MBlowerPI <- MBlower
  cover <- coverPI <- rep(NA,nsim)
  
  for(i in 1:nsim){
    datlist <- lapply(delta, function(d) simFunc(N, delta = d, psi = psi, zmax = zmax, zmin = zmin))
    
    PI1[i,] = unlist(lapply(datlist, function(d) mean(d$true.ymean.plus+(rnorm(N,e.mean1,1)/B) - d$true.ymean.min+(rnorm(N,e.mean2,1)/B))/
                              mean(d$true.amean.plus+(rnorm(N,e.mean1,1)/B) - d$true.amean.min+(rnorm(N,e.mean2,1)/B))))
    IF1[i,] = unlist(lapply(datlist, function(d) mean(f.num2(d))/mean(f.den2(d))))
    
    ll2 <- ul2 <- ul1 <- ll1 <- ll2PI <- ul2PI <- ll1PI <- ul1PI <- rep(NA,length(delta))
    for(jj in 1:length(delta)){
      psihat = IF1[i,jj]
      ifvals = (f.num2(datlist[[jj]]) - psihat*f.den2(datlist[[jj]])) / mean(datlist[[jj]]$true.amean.plus - datlist[[jj]]$true.amean.min)
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
                    emp.if = emp.cover[(length(emp.cover)/2+1):length(emp.cover)],
                    emp.pi = emp.cover[1:(length(emp.cover)/2)],
                    delta = rep(delta, length(K)),
                    rate = rep(round(K), each = length(delta)))
#write.csv(cover, file = 'simulationCoverageDub.csv')
#write.csv(df, file = 'simulationOutputDub.csv') 

library(ggplot2)
g<-ggplot(df) + 
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
#ggsave('type_rate_shift_plotsDub.png',g,height = 4, width = 7) #mu=2

g2<-ggplot(df) + 
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
#ggsave('type_rate_shift_plots_empSDDub.png',g2,height = 4, width = 7)

