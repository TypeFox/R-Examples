plot.nodeFillingOutput <-
function(x,...){
  if(class(x)!="nodeFillingOutput"){
       stop("Should be the output from the function nodeFilling")}
  nodes <- x$nodes
  poly <- x$poly
  hole.list <- x$hole.list
  plot(rbind(poly,poly[1,]),lwd=2,type="l",...)
  points(nodes,pch=19,cex=0.5,...)
  if(!is.null(hole.list)){
    number.holes <- length(hole.list)
    for(k in 1: number.holes){
      lines(rbind(hole.list[[k]],hole.list[[k]][1,]),lwd=2,...)
    }
  }
}
