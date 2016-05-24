# 21-09-2013: * Catch the input error, if there is only one row of probabilities in the triple case

plot.estPI <- function(x,col="black",highlight=NULL,hlCol="red",pch=20,zoom=FALSE,...){
  if(x$type=="single")
  {
    estPlotSingle(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  } else if(x$type=="pair")
  {
    estPlotPair(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  } else if(x$type=="triple")
  {
    if(!is.vector(x$obs)){
      if(nrow(x$probs)==1){
         stop("There is only one row of probabilities. You might want to use the option 'order=FALSE' in the estPI call in order to get all combinations.")
      }
    }
    estPlotTriple(x,col=col,highlight=highlight,hlCol=hlCol,pch=pch,zoom=zoom,...)
  }
}