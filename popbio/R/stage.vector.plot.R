stage.vector.plot<-function(stage.vectors, proportions=TRUE, legend.coords="topright", ylim=NULL, xlab="Years", ylab=NULL, col=rainbow(8), ... )
{
    p<-stage.vectors
    n<-dim(p)[1]                      #number of stage vectors
    if(is.null(n)){stop("stage.vectors should be a matrix with two or more stages")}
    x<-colnames(p)                      # x-axis names
    if(is.null(x)){x<-0:(dim(p)[2]-1)}  # start at 0,1,2,... if colnames missing
    if(length(col)<n){col<-rep(col,n)}  ## line colors (repeat if necessary)    
    if(proportions)                 ## plot proportions      
    {
       if(is.null(ylab)){ylab<- "Proportion in stage class" }
       p<-prop.table(p, 2)         ## Change counts to proportions using prop.table
       if(is.null(ylim)){ylim=c(min(p, na.rm=TRUE), max(p, na.rm=TRUE))}
       plot(x, p[1,], type='n', ylim=ylim, xlab=xlab, ylab=ylab,...)                
    }
    else                            ## OR plot total number 
    {
       if(is.null(ylab)){ylab <- "Number in stage class" }
       if(is.null(ylim)){ylim=c(floor(min(p, na.rm=TRUE)), ceiling(max(p, na.rm=TRUE))) }
       plot(x, p[1,], type='n', ylim=ylim, xlab=xlab, ylab=ylab, ... )      
    }
    ## order legend by decreasing order of mean number in stage vector
    y<-sort(apply(p,1,mean, na.rm=TRUE), index.return=TRUE, decreasing=TRUE)
    for (i in y$ix )                ## Loop through stage classes
    {                 
       lines(x, p[i,],lty=i, col=col[i], lwd=2)
    }
    leg.names<-paste(names(y$x), "")  ## pad names with trailing space for extra room
    if(leg.names[1]==" "){leg.names<-paste("row", y$ix, "")}
    legend(legend.coords[1],legend.coords[2], leg.names, lty=y$ix, col=col[y$ix], lwd=2)
}

