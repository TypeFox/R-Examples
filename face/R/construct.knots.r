construct.knots <- function(argvals,knots,knots.option,p){
  
if(length(knots)==1){
  allknots <- select.knots(argvals,knots,option=knots.option)
}

if(length(knots)>1){
  K = length(knots)-1 
  knots_left <- 2*knots[1]-knots[p:1+1]
  knots_right <- 2*knots[K] - knots[K-(1:p)]
  allknots <- c(knots_left,knots,knots_right)
}

return(allknots)

}