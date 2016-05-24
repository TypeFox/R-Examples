##################################################################################################
wrm.filter <- function (y, width, weight.type = 1, del = floor(width/2), extrapolate = TRUE) {
    if (is.numeric(y) == FALSE)
        stop("y must be numeric")
    len <- length(y)
    if (width > len) {
        cat("ERROR: width cannot exceed the length of the time series \n")
        width <- len
    }
    if (del >= width) {
        cat("ERROR: delay cannot exceed the window width \n")
        del <- width
    }
    

    if( width==0 ) { stop("invalid specification of 'width': value must be positive") } 
    
    if (weight.type!=0 & weight.type!=1 & weight.type!=2) { stop("invalid specification of 'weigth': possible values are 0, 1, 2 ") } 


    ## save the name of the input time series
    ts.name <- deparse(substitute(y))
    mis <- 0
    level <- rep(0, len)
    slope <- level
    n <- width
    xdat <- (del - n + 1):del
    we <- rep(1, width)
    if (weight.type == 1) {
        if(del == 0){
           we <- 1:n
        } else {
           we[1:(n-del)]   <- 1:(n-del)
           we[(n-del+1):n] <- n-del-(1:del)
        }     
    }
    if (weight.type == 2) {
    
        #Epanechnikov-Kern:
        epa <- function(xi,x = 0,h = 1){
        ifelse(abs(xi-x)<=h,1/h*3/4*(1- ((xi-x)/h)^2),0)
        }

        we=epa(xdat, x = 0, h = max(del, (n - del)))
    }
    for (x0 in (n - del):(len - del)) {
        if (is.na(y[x0 + del])) {
            y[x0 + del] <- level[x0 - 1] + slope[x0 - 1] * (del +1)
            mis <- 1
        }
        ydat <- y[x0 + xdat]
        erg <- WRMfit(xdat, ydat, 0, weight.vec = we)
        slope[x0] <- erg[2]
        level[x0] <- erg[1]
    }
    if (extrapolate == TRUE) {
        if ((n - del) > 1) {
            slope[1:(n - del - 1)] <- rep(slope[n - del], (n -
                del - 1))
            level[1:(n - del - 1)] <- level[n - del] + slope[n -
                del] * ((del - n + 1):(-1))
        }
        if (del > 0) {
            slope[(len - del + 1):len] <- rep(slope[len - del],
                del)
            level[(len - del + 1):len] <- level[len - del] +
                slope[len - del] * (1:del)
        }
    }
    else {
        if ((n - del) > 1) {
            slope[1:(n - del - 1)] <- NA
            level[1:(n - del - 1)] <- NA
        }
        if (del > 0) {
            slope[(len - del + 1):len] <- NA
            level[(len - del + 1):len] <- NA
        }
    }
    if (mis == 1) {
        cat("WARNING: Series contains missings \n")
    }
    weight.names <- c("Equal", "Triangular", "Epanechnikov")
    return(structure(list(y = y, level = level, slope = slope,
        del = del, width = width, weight.type = weight.names[weight.type+1], ts.name = ts.name),
        class = "wrm.filter"))
}



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default plot

plot.wrm.filter <- function(x, ...) {
    # Length of the time series
    N <- length(x$y)

    # Setting the y-limits
    ylims <- c(min(x$y,min(x$level,na.rm=TRUE),na.rm=TRUE),max(x$y,max(x$level,na.rm=TRUE),na.rm=TRUE))
    xlims <- c(1,N)
    
    # Defining the title
    t1 <- "Weighted Repeated Median Filter"
    titel <- ifelse(x$del==0,paste("Online ",t1,sep=""),t1)
    
    # Plot
    par(mar=c(4,4,4,7),oma=rep(0,4),mgp=c(2.5,1,0))
    plot(x$y,xlim=xlims,ylim=ylims,type="l",xlab="Time", ylab=x$ts.name, main=titel)
    lines(x$level, col="red" ,lwd=2)
    
    # Legend
    par(xpd=TRUE,cex=0.8)
    legend(par("usr")[2],mean(c(par("usr")[3],par("usr")[4])),c("Time Series","WRM"),xjust=0,yjust=0.5,lty=rep(1,2),lwd=c(1,2),col=c("black","red"),bty="n")
    par(xpd=FALSE,cex=1)
    }
    
    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Default output

print.wrm.filter <- function(x, ...) {
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
      cat(N-5," observations omitted \n")
    }
}
