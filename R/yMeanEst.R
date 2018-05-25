y.mean.est <- function(y,dat,algorithm,sl.lib = c("SL.gam","SL.glm","SL.glm.interaction","SL.glm","SL.ranger","SL.mean")){
  data = as.data.frame(cbind(y,dat))

  if(tolower(algorithm) == 'glm'){
    out = glm(y~., data = data, family = 'gaussian')
    }

  else if(tolower(algorithm) == 'superlearner'){
    library(SuperLearner)
    out = SuperLearner(Y=y, X=dat, SL.library=sl.lib, family=gaussian())
    }

  else if(tolower(algorithm) == 'ranger'){
    library(ranger)
    out = ranger::ranger(y~.,data = data, write.forest = T)
  }

  else if(tolower(algorithm) == 'random forest'){
    require(randomForest)
    out = randomForest(y~., data = data)
  }

  else{stop('Use ranger, superlearner or glm as algorithm')}

  return(out)
}
