plotRects <- function(R, grid=FALSE, grid.lty=3, grid.col="lightgray",
             density=NULL, angle=45, col=NA, 
             border=rainbow(nrow(R)), lty=1, lwd=1, add=FALSE, 
             xlim=NULL, ylim=NULL, xlab="", ylab="",
             main="", sub="")  
{
   if (!is.matrix(R) | ncol(R)!=4){
      stop("invalid argument 'R' - it must be a matrix with 4 columns")
   }
   if (sum(R[,1]>R[,2]) + sum(R[,3]>R[,4]) >=1){
      stop("invalid argument 'R' - R[i,1]<=R[i,2] or R[i,3]<=R[i,4] is 
        violated")
   }
   if (!is.logical(grid)){
      stop("invalid argument 'grid' - it must be logical")
   }
   if (!is.logical(add)){
      stop("invalid argument 'add' - it must be logical")
   }
  
   if (is.null(xlim)){
      xlim <- range(R[,1:2])
   }
   if (is.null(ylim)){
      ylim <- range(R[,3:4])
   }  

   if (!add)
   {
      # set up plotting region
      plot(c(0,0),type="n",xlim=xlim, ylim=ylim,xlab=xlab,ylab=ylab,
           main=main,sub=sub)  
   }

   # if grid==TRUE, plot grid lines
   if (grid){
      n <- nrow(R)
      x <- sort(R[,1:2])
      y <- sort(R[,3:4])
      abline(v=x[1:(2*n)], lty=grid.lty, col=grid.col)
      abline(h=y[1:(2*n)], lty=grid.lty, col=grid.col)    
   }

   # plot rectangles using the function rect()
   rect(R[,1],R[,3],R[,2],R[,4], density=density, angle=angle, col=col,
        border=border, lty=lty, lwd=lwd) 
}  

