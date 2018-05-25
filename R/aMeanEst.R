a.mean.est <- function(a,dat,algorithm,sl.lib = c("SL.gam","SL.glm","SL.glm.interaction","SL.glm","SL.ranger","SL.mean")){
  data = as.data.frame(cbind(a,dat))

  if(tolower(algorithm) == 'glm'){
    out = glm(a~., data = data, family = 'binomial')
    }

  else if(tolower(algorithm) == 'superlearner'){
    require(SuperLearner)
    out = SuperLearner(Y=a, X=dat, SL.library=sl.lib, family=binomial())}

  else if(tolower(algorithm) == 'ranger'){
    require(ranger)
    out = ranger::ranger(a~.,data = data, write.forest = T)
  }

  else if(tolower(algorithm) == 'random forest'){
    require(randomForest)
    out = randomForest(a~., data = data)
  }

  else{stop('Use ranger, superlearner, random forest or glm as algorithm')}

  return(out)
}
