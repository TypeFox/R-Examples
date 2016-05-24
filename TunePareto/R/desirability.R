
#desirability <- function(d1=0.01,x1,d2=0.99,x2)
#{
#  b1 <- (-log(-log(d2))+log(-log(d1)))/(x2-x1)
#  b0 <- -log(-log(d1))-b1*x1
#  return(c(b0,b1))
#}

#calculateDesirabilities <- function(objectiveValues, desirabilities)
#{
#  d <- sapply(1:ncol(objectiveValues),function(objective)
#       {
#        vals <- objectiveValues[,objective]
#        des <- desirabilities[[objective]]
#        return(exp(-exp(-(des[1]+des[2]*vals))))
#       })
# return(apply(d,1,function(x)prod(x)^(1/length(x))))
#}

rankByDesirability <- function(tuneParetoResult, desirabilityIndex, optimalOnly=TRUE)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")


  if (optimalOnly)
    objectiveValues <- tuneParetoResult$bestObjectiveValues
  else
    objectiveValues <- tuneParetoResult$testedObjectiveValues
    
  desirabilities <- apply(objectiveValues,1,desirabilityIndex)
  idx <- order(desirabilities, decreasing=TRUE)
  
  mat <- data.frame(objectiveValues,desirabilities)
  colnames(mat)[ncol(mat)] <- "Desirability"
  
  mat <- mat[idx,]
  return(mat)
}
