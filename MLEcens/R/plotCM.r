plotCM <- function(cm, col=c("white","black"), at.x=NULL, at.y=NULL, 
       xlab="Observation rectangles", ylab="Maximal intersections", 
       main="", sub=""){
 
   if (!is.matrix(cm)){
       stop("invalid argument 'cm'")
   }  
 
   image(c(1:ncol(cm)),c(1:nrow(cm)),t(cm),col=col,xlab=xlab,ylab=ylab, 
         main=main,sub=sub,axes=FALSE) 
   axis(1, at=at.x)
   axis(2, at=at.y)
   box()
}

