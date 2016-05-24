mogavs.formula <-
function(formula, data, maxGenerations=10*ncol(x),popSize=ncol(x),noOfOffspring=ncol(x),crossoverProbability=0.9,mutationProbability=1/ncol(x),kBest=1,plots=F,additionalPlots=F, ...){

  mf<-model.frame(formula=formula, data=data)
  x<- model.matrix(attr(mf, "terms"), data=mf)[,-1] #remove the intercept
  y <- model.response(mf)
  
  est<-mogavs.default(x, y, maxGenerations,popSize,noOfOffspring,crossoverProbability,mutationProbability,kBest,plots, additionalPlots)
  
  est$call <- match.call()
  est$formula <- formula
  est
}
