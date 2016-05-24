#'@export
plot.DStree<-function(x,prob="surv",select=NULL,...){
  
  if (!prob %in% c("haz","surv")){
    stop("argument 'prob' must be either 'haz' or 'surv'")
  }
  
  fitr<-x
  fitr$method<-"anova"
  class(fitr)<-"rpart"
  
  n.lev <- length(fitr$parms[[1]]) 
  
  maxsurv <- x$parms[[1]][length(x$parms[[1]])]
  med<-x$frame$yval
  med[which(is.na(med))] <- paste(">",maxsurv)
  node.fun1 <- function(x, labs, digits, varlen)
  {
    paste("Leaf: ",rownames(x$frame),"\nSurv. Median:", med,
          "\n\nn=",x$frame$n)
  }
  
  
  ind <- which(fitr$frame$var=="<leaf>")
  if (missing(select)){
    ind.select <- ind
  }else
  {
    ind.select <- which(rownames(fitr$frame) %in% select)
  }
  
  if (length(ind.select)==0){
    stop("selected leaves in argument 'select' are not part of the fit")
  }
  
  ind.leaf<-intersect(ind.select,ind)
  
  n.leaf<-length(ind.leaf)
  
  if(prob=="surv"){
    cols <- (1+n.lev):(2*n.lev)
    main <- "Survival Probabilities"
  }
  else{
    cols <- 1:(n.lev)
    main <- "Hazard Probabilities"
  }
  
  par(xpd=TRUE)
  plot(fitr$frame$yval2[ind.leaf[1],cols],type="S",xlab="Time",
       ylab="", main=main, col= rainbow(n.leaf)[1],ylim=c(0,1))
  
  for(i in 2:n.leaf){
    lines(fitr$frame$yval2[ind.leaf[i],cols],type="S",col= rainbow(n.leaf)[i])
  }
  legend("topright", legend = rownames(x$frame)[ind.leaf], 
         col= rainbow(n.leaf)[1:n.leaf], lty = 1, pt.cex = 1,
         cex=0.5,pch = 1.5)
  
  if(missing(select)){
  prp(fitr,node.fun=node.fun1,...)}
}