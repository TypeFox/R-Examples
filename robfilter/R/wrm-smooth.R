
wrm.smooth <- function(x,y, h, xgrid,  weight=2){
  if (weight!=1 & weight!=2 & weight!=3 & weight!=4 & weight!=5) { stop("invalid specification of 'weigth': possible values are 1, 2, 3, 4, 5 ") } 

  if (sum(is.na(cbind(x,y)))){
      new <-  na.omit(cbind(x,y))
      x <- new[,1]
      y <- new[,2]
      cat("WARNING: Missing values omitted. \n")
    }
  if (missing(h)){
    stop("No automatic bandwidth selection routine implemented yet. Please specify bandwidth manually. \n") 
  }
  if (missing(xgrid)){
     if (length(x)<=100) xgrid <- sort(x) else xgrid <- seq(min(x),max(x), l=100)
  }
  N     <- length(xgrid)
  level <- slope <- rep(0,N)
  for (i in 1:N){
     weights   <- kernels(x, xgrid[i], h, weight)
     fitted    <- WRMfit(x,y, xgrid[i], weights)
     level[i]  <- fitted[1]
     slope[i]  <- fitted[2]
  } 
  weight.names <- c("Triangular", "Epanechnikov", "Gaussian", "Biweight", "Uniform")

  return(structure(list(x=x, 
        y = y,
            level=level, 
            slope=slope, 
            h=h, 
            xgrid=xgrid,
            weight=weight.names[weight]), class="wrm.smooth"))
}





# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot

plot.wrm.smooth <- function(x, ...) {
    # Length of the time series
    N <- length(x$y)
    # Weight names
    weight <- x$weight
    
    # Setting the y-limits
    ylims <- c(min(x$y,min(x$level,na.rm=TRUE),na.rm=TRUE),max(x$y,max(x$level,na.rm=TRUE),na.rm=TRUE))
    xlims <- c(min(x$x),max(x$x))
    
    # Defining the title
    titel <- "Weighted Repeated Median Smoother"
    
    # Plot
    par(mar=c(4,4,4,8),oma=rep(0,4),mgp=c(2.5,1,0))
    plot(x$x,x$y, type="p", xlim=xlims,ylim=ylims,xlab="x", ylab="y", main=titel)
    lines(x$xgrid,x$level, col="red" ,lwd=2)
    
    # Legend
    par(xpd=TRUE,cex=0.8)
    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c("Original","WRM Smoother"),xjust=0,yjust=0.5,lty=rep(1,2),lwd=c(1,2),col=c("black","red"),bty="n")
    par(xpd=FALSE,cex=1)
    }


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output
print.wrm.smooth <- function(x, ...) {
    # length of the input time series
    N <- length(x$y)
    if(N <= 50){
      print(list(level=round(x$level,6), slope=round(x$slope,6), sigma=round(x$sigma,6)))
    } else {
      cat("$level")
      L <- data.frame("[1]",t(round(x$level[1:5],6)),"...")
      dimnames(L)[[1]] <- c(" ")
      dimnames(L)[[2]] <- c(" ","  ","   ", "    ","     ","      ","       ")
      print(L[1,])
      cat("$slope")
      Sl <- data.frame("[1]",t(round(x$slope[1:5],6)),"...")
      dimnames(Sl)[[1]] <- c(" ")
      dimnames(Sl)[[2]] <- c(" ","  ","   ", "    ","     ","      ","       ")
      print(Sl[1,])
    }
}
