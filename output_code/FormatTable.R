dubshift <- read.csv("dubshiftEven_5foldRR.csv")[,-1]
dubshift.child <- read.csv("dubshiftEvenChild_5foldRR.csv")[,-1]
to.print <- rbind(delta = round(dubshift[7,]),
                  all.visits = round(dubshift[1,],3),
                  unif.ci = apply(dubshift[5:6,],2, function(x) paste("(",paste(round(x,3), collapse = ","),")", sep = "")),
                  compliers = round(dubshift[9,],3),
                  child.visits = round(dubshift.child[1,],3),
                  unif.ci.child = apply(dubshift.child[5:6,],2, function(x) paste("(",paste(round(x,3), collapse = ","),")", sep = "")),
                  compliers.child = round(dubshift.child[9,],3)
                  )
library(xtable)
xtable(to.print)
