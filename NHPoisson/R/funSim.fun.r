funSim.fun <-
function(i,lambda,fun.name,fun.args=NULL)
{
posNH<-simNHP.fun(lambda=lambda)$posNH
  fun.args <- c(list(posNH), fun.args)
aux<-do.call(fun.name,fun.args)
return(aux)
}
