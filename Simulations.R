#################################################################################
#     this file includes all code to produce simulations                        #
#     use this to reproduce figures 1 and 2 in the paper                        #
#################################################################################

############## positivity simulation (figure 1) ###########
rm(list = ls())

today = format(Sys.time(), "%Y%m%d")
set.seed(123)
library(ggplot2)
library(truncnorm)

# generate covariate and postivity limits
N <- 5000
x <- rbinom(N, 1, .5)
mu <- 2*x - 1
Upp.lim <- 1*(x==0) + 3*(x==1)
Low.lim <- -3*(x==0) + -1*(x==1)
Upp.lim0 <- Upp.lim[x==0][1]; Low.lim0 <- Low.lim[x==0][1]
Upp.lim1 <- Upp.lim[x==1][1]; Low.lim1 <- Low.lim[x==1][1]

# generate truncated normal
z <- rep(-5, N)
for(i in 1:N){while(  (z[i]<Low.lim[i]) | (z[i]>Upp.lim[i])  ){z[i] <- rnorm(1, mu[i], .5)}}

delta = .1
g2<-ggplot()+
  geom_histogram(aes(z[x==1], fill = 'X=1'), alpha = .7, bins = 50)+
  geom_histogram(aes(z[x==0], fill = 'X=0'), alpha = .7, bins = 50)+
  theme_bw()+
  geom_rect(aes(xmin=Upp.lim0-delta, xmax=Upp.lim0, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.7)+
  geom_rect(aes(xmin=Low.lim0, xmax=Low.lim0+delta, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.7)+
  geom_rect(aes(xmin=Low.lim1, xmax=Low.lim1+delta, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.7)+
  geom_rect(aes(xmin=Upp.lim1-delta, xmax=Upp.lim1, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.7)+
  geom_rect(aes(xmin=Upp.lim0, xmax=Upp.lim1, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.2)+
  geom_rect(aes(xmin=Low.lim0, xmax=Low.lim1, ymin=-1, ymax=length(z)/100), color = 'black',alpha=0.2)+
  xlab("Z") + ylab("Frequency")+
  scale_fill_grey(start = 0, end = .5)+
  guides(fill=FALSE)+
  ggtitle("Low Overlap")
ggsave(filename = 'PosOverlap4bars.png', g2, width = 7, height = 4)

# check how likely positivity violations are
length(which(z[x==0]+delta>Upp.lim0))+length(which(z[x==0]-delta<Low.lim0))+
  length(which(z[x==1]+delta>Upp.lim1))+length(which(z[x==1]-delta<Low.lim1))
1-ptruncnorm(Upp.lim0-delta,Low.lim0,Upp.lim0,min(mu),.5)


############## rate simulations (figure 2) ###########
rm(list = ls())

today = format(Sys.time(), "%Y%m%d")
set.seed(123)
library(ggplot2)
library(AER)
library(gridExtra)

# errors to add to nuisance parameters -- change these to test the effects of different error terms
mu.error <- function(z){rnorm(N, 2*z + 1, 1)/B}
la.error <- function(z){rnorm(N, 2*z + 1, 1)/B}
rat.error <- function(z){rnorm(N,1,1)/B}
n.sim = 500
deltas = seq(.5,4,length.out = 15)
K = c(1.99,2.99,3.99,5.99)
Ns = c(100,1000,5000,10000)


# functions to simulate and estimate
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
f.num <- function(df){
  mu0 = df$true.ymean + mu.error(df$z)
  muP = df$true.ymean.plus + mu.error(df$z+delta)
  muM = df$true.ymean.min + mu.error(df$z-delta)
  ratM = df$true.z.min/df$true.z + rat.error(df$z-delta)
  ratP = df$true.z.plus/df$true.z + rat.error(df$z+delta)
  (ratM*(df$y - mu0) + muP) - (ratP*(df$y - mu0) + muM)
}
f.den <- function(df){
  la0 = df$true.amean + la.error(df$z)
  laP = df$true.amean.plus + la.error(df$z+delta)
  laM = df$true.amean.min + la.error(df$z-delta)
  ratM = df$true.z.min/df$true.z + rat.error(df$z-delta)
  ratP = df$true.z.plus/df$true.z + rat.error(df$z+delta)
  (ratM*(df$a - la0) + laP) - (ratP*(df$a - la0) + laM)
}
pi.num <- function(df){(df$true.ymean.plus+mu.error(df$z+delta)) - (df$true.ymean.min + mu.error(df$z-delta))}
pi.den <- function(df){(df$true.amean.plus + la.error(df$z+delta) ) - (df$true.amean.min + la.error(df$z-delta))}
get_if <- function(df){mean(f.num(df))/mean(f.den(df))}
get_pi <- function(df){mean(pi.num(df))/mean(pi.den(df))}
get_tsls_est <- function(df){summary(ivreg(y~a+X1+X2+X3+X3 | z+X1+X2+X3+X3, data = df))$coef["a","Estimate"]}
get_tsls_sd <- function(df){summary(ivreg(y~a+X1+X2+X3+X3 | z+X1+X2+X3+X3, data = df))$coef["a","Std. Error"]}
psi_sd <- function(df){sd((f.num(df) - get_if(df)*f.den(df))/mean((df$true.amean.plus+la.error(df$z+delta))-(df$true.amean.min+la.error(df$z-delta))))/sqrt(dim(df)[1])}
single.delta <- function(delta){
  dfs = lapply(1:n.sim, function(x) simFunc(N=N, delta = delta))
  if.ests = unlist(lapply(dfs, function(df) get_if(df)))
  pi.ests = unlist(lapply(dfs, function(df) get_pi(df)))
  th.sds = unlist(lapply(dfs, function(df) psi_sd(df)))
  tsls.ests = unlist(lapply(dfs, function(df) get_tsls_est(df)))
  tsls.sds = unlist(lapply(dfs, function(df) get_tsls_sd(df)))
  return(data.frame(if.ests = if.ests, pi.ests = pi.ests, th.sds = th.sds, tsls.ests = tsls.ests, tsls.sds = tsls.sds))
}

# run the simulation across all sample sizes and shift values
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

# save the output
#lapply(output, function(k) write.csv(k,file=paste("~/Dropbox/double robust causality/df_",k$s.size[1],"_",today,".csv",sep="")))

# make figures
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
plots <- lapply(output, function(x) make.my.plot(x[x$type != "TSLS",]))
#for(i in 1:length(plots)){ggsave(plot = plots[[i]], filename = paste("~/Dropbox/double robust causality/Figures/simulation_plot_",Ns[i],today,".png", sep = ""), height = 4, width = 7)}
all.output <- do.call(rbind.data.frame, output)
combo.plot <- ggplot(all.output[(all.output$type!="TSLS")&(all.output$s.size!=5000),]) +
  geom_hline(yintercept = psi, colour = 'red')+
  geom_point(aes(x = delta.val, y = mean, shape = type))+
  geom_line(aes(x = delta.val, y = mean, colour = type))+
  geom_ribbon(aes(ymin = lower.emp, ymax = upper.emp, x = delta.val, fill = type), alpha = .3) +
  facet_wrap(~s.size + rate.val)+
  theme_bw() +
  scale_color_manual(values = c('black', 'grey'))+
  scale_fill_manual(values = c('black', 'grey'))+
  theme(legend.position="bottom")+
  coord_cartesian(ylim = c(0, 4)) +
  labs(y = paste('Estimates (',n.sim,' simulations)', sep = ""), x = "Shift Amount",
       title = "Estimates by estimator type, error rate and shift amount")
#ggsave(plot = combo.plot, filename = paste("~/Dropbox/double robust causality/Figures/simulation_plot_combo_",today,".png", sep = ""), height = 8, width = 7)


make.my.plot.cf <- function(df){
  # use closed form CI's
  ggplot(df) +
    geom_hline(yintercept = psi, colour = 'red')+
    geom_point(aes(x = delta.val, y = mean, shape = type))+
    geom_line(aes(x = delta.val, y = mean, colour = type))+
    geom_ribbon(aes(ymin = lower.th, ymax = upper.th, x = delta.val, fill = type), alpha = .3) +
    facet_wrap(~rate.val)+
    theme_bw() +
    scale_color_manual(values = c('black', 'grey'))+
    scale_fill_manual(values = c('black', 'grey'))+
    theme(legend.position="bottom")+
    coord_cartesian(ylim = c(0, 4)) +
    labs(y = paste('Estimates (',n.sim,' simulations, sample size = ', df$s.size,')', sep = ""), x = "Shift Amount",
         title = "Estimates by estimator type, error rate and shift amount (closed form SD)")
}
plots.cf <- lapply(output, function(x) make.my.plot.cf(x[x$type != "TSLS",]))
#for(i in 1:length(plots)){ggsave(plot = plots.cf[[i]], filename = paste("~/Dropbox/double robust causality/Figures/simulation_plot_",Ns[i],today,"_cf.png", sep = ""), height = 4, width = 7)}
