
plot_ldistr <- function(x,...)
{


layout(matrix(1:length(x$QUAD),nrow=1))

for(i in 1:length(x$QUAD))
{
  
plot(x$QUAD[[i]]$nodes,x$QUAD[[i]]$weights, type="h",xlab="Quadrature points", ylab="", main=paste("Group",names(x$ZLpar)[i]), ...)

}


}










