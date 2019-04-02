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
covs0 = cbind(white,loslastloc,
             county.dummies,
             ageyrs,
             urban,
             priorarrests,married,violent,
             lsirscore,
             highschoolgrad,
             custody_level,
             numofpriorinc,
             mh,dat[,c(startdummy:(enddummy-1))],
             numoftotalmisconducts)

# using minimum distance prison as a dummy
covs = cbind(white,loslastloc,
             #county.dummies,
             minName.dummies,
             ageyrs,
             urban,
             priorarrests,married,violent,
             lsirscore,
             highschoolgrad,
             custody_level,
             numofpriorinc,
             mh,dat[,c(startdummy:(enddummy-1))],
             numoftotalmisconducts)

# using minimum distance prison (no dummies)
covs3 = cbind(white,loslastloc,
             #county.dummies,
             minTimeName,
             ageyrs,
             urban,
             priorarrests,married,violent,
             lsirscore,
             highschoolgrad,
             custody_level,
             numofpriorinc,
             mh,dat[,c(startdummy:(enddummy-1))],
             numoftotalmisconducts)

source('~jacquelinemauro/Dropbox/double robust causality/tempR/longlat.R')
# using minimum distance prison longitude and latitude 
min.prison <- apply(dat[,2:26],1,which.min)
min.lat <- locs[min.prison,1]
min.long <- locs[min.prison,2]
covs4 = cbind(white,loslastloc,
              min.lat,min.long,
              ageyrs,
              urban,
              priorarrests,married,violent,
              lsirscore,
              highschoolgrad,
              custody_level,
              numofpriorinc,
              mh,dat[,c(startdummy:(enddummy-1))],
              numoftotalmisconducts)

detach(dat)
names(covs0) <- sapply(c(1:dim(covs0)[2]), function(k) paste('x',k,sep = ""))
names(covs) <- sapply(c(1:dim(covs)[2]), function(k) paste('x',k,sep = ""))
names(covs3) <- sapply(c(1:dim(covs3)[2]), function(k) paste('x',k,sep = ""))
names(covs4) <- sapply(c(1:dim(covs4)[2]), function(k) paste('x',k,sep = ""))

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

#### run double shift est  ####
delta1 = c(20,30,60,90)
delta2 = c(0,20,30,60)
output <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'flexcode',
                             zmin = 0,zmax = 500, nfolds = 21, pos.cutoff = 100)
p <- plot.cace(output)
ggsave(filename = 'recidivismEs_1fold.png',p,height = 4, width = 7)

delta = c(20,30,60,90)
output <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output) + 
  geom_hline(yintercept = 0) + 
  theme_bw()+
  ggtitle('Increase in recidivism when visitation is lost')
ggsave(filename = 'recidivismEs_2foldNoDummsRangRangGLM.png',p,height = 4, width = 7)

output.child <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.childvisits,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'flexcode',
                             zmin = 0,zmax = 500, nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output.child)
ggsave(filename = 'recidivismEst_child.png',p,height = 4, width = 7)

output.child <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                             a = dat$no.childvisits,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output.child) + 
  geom_hline(yintercept = 0) + 
  theme_bw()+
  ggtitle('Increase in recidivism when child visitation is lost')
ggsave(filename = 'recidivismEs_child2foldNoDummsRangRangGLM.png',p,height = 4, width = 7)

#### run single shift est, shifting up ####
delta = c(20,40,60,90,120)
# one fold reduces variance compared to 2
output <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 1, pos.cutoff = Inf)
p1 <- plot.cace(output) + 
  ggtitle("Effect of No Visitation on Recidivism, Moving Away")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftUp.png',p1,height = 4, width = 7)

out1df<- t(matrix(unlist(output[-5]), ncol = length(delta), byrow = T))
write.csv(out1df, 'recidivismEst_SingleShiftUp.csv')

# with superlearner for means and ranger in the kernel
output <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 1, pos.cutoff = Inf)
p1 <- plot.cace(output) + 
  ggtitle("Effect of No Visitation on Recidivism, Moving Away")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftUpSL.png',p1,height = 4, width = 7)

out1df<- t(matrix(unlist(output[-5]), ncol = length(delta), byrow = T))
write.csv(out1df, 'recidivismEst_SingleShiftUpSL.csv')

# with superlearner for means and ranger in the kernel, 90% unifband
output <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 1, pos.cutoff = Inf, alpha = 0.1)
p1 <- plot.cace(output) + 
  ggtitle("Effect of No Visitation on Recidivism, Moving Away")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftUpSL90.png',p1,height = 4, width = 7)

out1df<- t(matrix(unlist(output[-5]), ncol = length(delta), byrow = T))
write.csv(out1df, 'recidivismEst_SingleShiftUpSL90.csv')

# with factors, not dummies
# one fold reduces variance compared to 2
output <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 1, pos.cutoff = Inf)
p1 <- plot.cace(output) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Away")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftUpNoDumms.png',p1,height = 4, width = 7)



output.child <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                   nfolds = 1, pos.cutoff = Inf)
p2 <- plot.cace(output.child) + 
  ggtitle("Effect of No Child Visitation on Recidivism, Moving Away")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftUp.png',p2,height = 4, width = 7)

out2df<- t(matrix(unlist(output.child[-5]), ncol = length(delta), byrow = T))
write.csv(out2df, 'recidivismEst_childSingleShiftUp.csv')

# with superlearner/ranger
output.child <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                   nfolds = 1, pos.cutoff = Inf)
p2 <- plot.cace(output.child) + 
  ggtitle("Effect of No Child Visitation on Recidivism, Moving Away")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftUpSL.png',p2,height = 4, width = 7)

out2df<- t(matrix(unlist(output.child[-5]), ncol = length(delta), byrow = T))
write.csv(out2df, 'recidivismEst_childSingleShiftUpSL.csv')

output.child <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                   nfolds = 1, pos.cutoff = Inf)
p2 <- plot.cace(output.child) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Away")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftUpNoDumms.png',p2,height = 4, width = 7)


#### run single shift est, shifting down ####
delta = -c(20,40,60,90,120)
output2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$visitslastlocyn1,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime),
                             nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Closer") +
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftDown.png',p3,height = 4, width = 7)

out3df<- t(matrix(unlist(output2[-5]), ncol = length(delta), byrow = T))
write.csv(out3df, 'recidivismEst_SingleShiftDown.csv')

output2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$visitslastlocyn1,delta = delta,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime),
                              nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Closer") +
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftDownSL.png',p3,height = 4, width = 7)

out3df<- t(matrix(unlist(output2[-5]), ncol = length(delta), byrow = T))
write.csv(out3df, 'recidivismEst_SingleShiftDownSL.csv')

output2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                              a = dat$visitslastlocyn1,delta = delta,
                              Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime),
                              nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Closer") +
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftDownNoDumms.png',p3,height = 4, width = 7)

output.child2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                   a = dat$childyn,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                   nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Closer")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftDown.png',p4,height = 4, width = 7)

out4df<- t(matrix(unlist(output.child2[-5]), ncol = length(delta), byrow = T))
write.csv(out4df, 'recidivismEst_childSingleShiftDown.csv')

output.child2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                    a = dat$childyn,delta = delta,
                                    Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                                    zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                    nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Closer")+
  theme_bw()+
  xlab("Shift amount (driving minutes)")+
  ylab("Change in Recidivism Probability")+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftDownSL.png',p4,height = 4, width = 7)

out4df<- t(matrix(unlist(output.child2[-5]), ncol = length(delta), byrow = T))
write.csv(out4df, 'recidivismEst_childSingleShiftDownSL.csv')

output.child2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                                    a = dat$childyn,delta = delta,
                                    Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                    zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                    nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Closer")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftDownNoDumms.png',p4,height = 4, width = 7)

#### run using Reinc as y ----
delta = c(20,40,60,90,120)
# one fold reduces variance compared to 2
outputReinc <- single.shift.range(y = dat$NCReinc3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 1, pos.cutoff = Inf)
p1 <- plot.cace(outputReinc) + 
  ggtitle("Effect of All Visitation on Re-Incarceration, Moving Away")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'reincEst_SingleShiftUp.png',p1,height = 4, width = 7)

write.csv(t(matrix(unlist(outputReinc[-5]), ncol = length(delta), byrow = T)), 'reincEst_SingleShiftUp.csv')


output.childReinc <- single.shift.range(y = dat$NCReinc3,z = dat$total_time,x = covs,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                   nfolds = 1, pos.cutoff = Inf)
p2 <- plot.cace(output.childReinc) + 
  ggtitle("Effect of Child Visitation on Re-Incarceration, Moving Away")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'reincEst_childSingleShiftUp.png',p2,height = 4, width = 7)

write.csv(t(matrix(unlist(output.childReinc[-5]), ncol = length(delta), byrow = T)), 'reincEst_childSingleShiftUp.csv')

delta = -c(20,40,60,90,120)
output2Reinc <- single.shift.range(y = dat$NCReinc3,z = dat$total_time,x = covs,
                              a = dat$visitslastlocyn1,delta = delta,
                              Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime),
                              nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2Reinc) + 
  ggtitle("Effect of All Visitation on Re-Incarceration, Moving Closer") +
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'reincEst_SingleShiftDown.png',p3,height = 4, width = 7)

write.csv(t(matrix(unlist(output2Reinc[-5]), ncol = length(delta), byrow = T)), 'reincEst_SingleShiftDown.csv')

output.child2Reinc <- single.shift.range(y = dat$NCReinc3,z = dat$total_time,x = covs,
                                    a = dat$childyn,delta = delta,
                                    Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                    zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                    nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2Reinc) + 
  ggtitle("Effect of Child Visitation on Re-Incarceration, Moving Closer")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'reincEst_childSingleShiftDown.png',p4,height = 4, width = 7)

write.csv(t(matrix(unlist(output.child2Reinc[-5]), ncol = length(delta), byrow = T)), 'reincEst_childSingleShiftDown.csv')

#### run single shift est, shifting down, a = no visits ####
delta = -c(20,40,60,90,120)
output2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.visits.last.loc,delta = delta,
                              Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime),
                              nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Closer") +
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftDown.png',p3,height = 4, width = 7)

output2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                              a = dat$visitslastlocyn1,delta = delta,
                              Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime),
                              nfolds = 1, pos.cutoff = Inf)
p3 <- plot.cace(output2) + 
  ggtitle("Effect of All Visitation on Recidivism, Moving Closer") +
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_SingleShiftDownNoDumms.png',p3,height = 4, width = 7)

output.child2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                    a = dat$childyn,delta = delta,
                                    Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                    zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                    nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Closer")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftDown.png',p4,height = 4, width = 7)

output.child2 <- single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                                    a = dat$childyn,delta = delta,
                                    Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                    zmin = min(dat$minTime),zmax = max(dat$maxtime),
                                    nfolds = 1, pos.cutoff = Inf)
p4 <- plot.cace(output.child2) + 
  ggtitle("Effect of Child Visitation on Recidivism, Moving Closer")+
  theme_bw()+
  geom_hline(yintercept = 0)
ggsave(filename = 'recidivismEst_childSingleShiftDownNoDumms.png',p4,height = 4, width = 7)

### put them all together ---
up.all <- read.csv('~jacquelinemauro/Dropbox/double robust causality/recidivismEst_SingleShiftUp.csv')[,-1]
up.chi <- read.csv('~jacquelinemauro/Dropbox/double robust causality/recidivismEst_childSingleShiftUp.csv')[,-1]
dn.all <- read.csv('~jacquelinemauro/Dropbox/double robust causality/recidivismEst_SingleShiftDown.csv')[,-1]
dn.chi <- read.csv('~jacquelinemauro/Dropbox/double robust causality/recidivismEst_childSingleShiftDown.csv')[,-1]

all <- rbind(up.all, up.chi, dn.all, dn.chi)
names(all) <- c('est','sd','llp','ulp','llu','ulu','delta','compliers')
all$direction <- c(rep('Farther',10), rep('Closer',10))
all$type <- c(rep(c(rep('All',5), rep('Child',5)),2))
all$delta <- abs(all$delta)
p <- ggplot(all) + geom_line(aes(x = delta, y = est, colour = "Estimate"))+
  geom_point(aes(x = delta, y = est, colour = "Estimate"))+
  geom_ribbon(aes(ymin = llp, ymax = ulp, x = delta, fill = "Pointwise"), alpha = .3) +
  geom_ribbon(aes(ymin = llu, ymax = ulu, x = delta, fill = "Uniform"), alpha = .3) +
  facet_grid(type~direction, scales = 'free')+
  theme_bw()+
  xlim(c(20,120))+ geom_hline(yintercept = 0)+
  labs(x = expression(delta), y = expression(psi), title = "Change in Recidivism among Compliers",
       fill = 'band type', colour = '')
ggsave(filename = 'recidivismEst_all.png',p,height = 4, width = 7)


up.all <- read.csv('recidivismEst_SingleShiftUpSL.csv')[,-1]
up.chi <- read.csv('recidivismEst_childSingleShiftUpSL.csv')[,-1]
dn.all <- read.csv('recidivismEst_SingleShiftDownSL.csv')[,-1]
dn.chi <- read.csv('recidivismEst_childSingleShiftDownSL.csv')[,-1]

all <- rbind(up.all, up.chi, dn.all, dn.chi)
names(all) <- c('est','sd','llp','ulp','llu','ulu','delta','compliers')
all$direction <- c(rep('farther',10), rep('closer',10))
all$type <- c(rep(c(rep('all',5), rep('child',5)),2))
all$delta <- abs(all$delta)
p <- ggplot(all) + geom_line(aes(x = delta, y = est, colour = "Estimate"))+
  geom_point(aes(x = delta, y = est, colour = "Estimate"))+
  geom_ribbon(aes(ymin = llp, ymax = ulp, x = delta, fill = "Pointwise"), alpha = .3) +
  geom_ribbon(aes(ymin = llu, ymax = ulu, x = delta, fill = "Uniform"), alpha = .3) +
  facet_grid(type~direction, scales = 'free')+
  theme_bw()+
  xlim(c(20,120))+ geom_hline(yintercept = 0)+
  labs(x = expression(delta), y = expression(psi), title = "Change in Recidivism among Compliers",
       fill = 'band type', colour = '')
ggsave(filename = 'recidivismEst_allSL.png',p,height = 4, width = 7)


#### run double shift est   ####
# check that we recover the single shift
# tried using min prison lat/long
# instead of matrix of 0/1's. this makes it run a lot faster. but CI's wider.
delta1 <- c(20,40,60,90,120)
delta2 <- rep(0,length(delta1))
delta = delta1
temp = single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                          a = dat$no.visits.last.loc,delta = delta,
                          Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                          zmin = min(dat$minTime),zmax = max(dat$maxtime),
                          nfolds = 1, pos.cutoff = Inf)
temp2 = single.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs4,
                          a = dat$no.visits.last.loc,delta = delta,
                          Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                          zmin = min(dat$minTime),zmax = max(dat$maxtime),
                          nfolds = 1, pos.cutoff = Inf)
write.csv(matrix(unlist(temp2[-5]),ncol = length(delta),byrow = T),file = 'deleteme.csv')
temp3 = double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs4,
                           a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                           Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                           zmin = min(dat$total_time),zmax = max(dat$total_time),
                           nfolds = 2, pos.cutoff = Inf)
write.csv(matrix(unlist(temp3[-5]),ncol = length(delta),byrow = T),file = 'deleteme2.csv')

# run for reals

# first, have delta1 = -delta2
delta1 <- c(20,40,60,90,120)
delta2 <- -delta1
output1 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                             a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                             Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 2, pos.cutoff = 100)
p <- plot.cace.double(output1)
write.csv(matrix(unlist(output1[-5]),ncol = length(delta1),byrow = T),file = 'dubshiftEven.csv')
ggsave(plot = p, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShiftEven.png',height = 4, width = 7)

# child visitation delta1 = -delta2 
delta1 <- c(20,40,60,90,120)
delta2 <- -delta1
output1c <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.childvisits,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 2, pos.cutoff = 100)
p <- plot.cace.double(output1c)
write.csv(matrix(unlist(output1c[-5]),ncol = length(delta1),byrow = T),file = 'dubshiftEvenChild.csv')
ggsave(plot = p, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShiftEvenChild.png',height = 4, width = 7)


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
delta1 <- c(-60,0,60,120)
delta2 <- c(-120,-60,0,60)
output4 <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                              a = dat$no.visits.last.loc,delta1 = delta1,delta2=delta2,
                              Y.est = 'superlearner', A.est = 'superlearner',Z.est = 'ranger',
                              zmin = min(dat$minTime),zmax = max(dat$maxtime), nfolds = 2, 
                              pos.cutoff = 50)
p4 <- plot.cace.double(output4)
write.csv(matrix(unlist(output4[-5]),ncol = length(delta1),byrow = T),file = 'dubshift60Window.csv')
ggsave(plot = p4, filename = '~jacquelinemauro/Dropbox/double robust causality/Figures/DubShift60Window.png',height = 4, width = 7)


# not done
delta = c(20,30,60,90)
output <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                             a = dat$no.visits.last.loc,delta = delta,
                             Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                             zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                             nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output) + 
  geom_hline(yintercept = 0) + 
  theme_bw()+
  ggtitle('Increase in recidivism when visitation is lost')
ggsave(filename = 'recidivismEs_2foldNoDummsRangRangGLM.png',p,height = 4, width = 7)

output.child <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'flexcode',
                                   zmin = 0,zmax = 500, nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output.child)
ggsave(filename = 'recidivismEst_child.png',p,height = 4, width = 7)

output.child <- double.shift.range(y = dat$NCRecid3,z = dat$total_time,x = covs3,
                                   a = dat$no.childvisits,delta = delta,
                                   Y.est = 'ranger', A.est = 'ranger',Z.est = 'glm',
                                   zmin = min(dat$minTime),zmax = max(dat$maxtime), 
                                   nfolds = 2, pos.cutoff = 100)
p <- plot.cace(output.child) + 
  geom_hline(yintercept = 0) + 
  theme_bw()+
  ggtitle('Increase in recidivism when child visitation is lost')
ggsave(filename = 'recidivismEs_child2foldNoDummsRangRangGLM.png',p,height = 4, width = 7)
