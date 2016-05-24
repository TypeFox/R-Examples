## .................................................................
summary.structureSim <- function(object, index=c(1:15), eigenSelect=NULL, ...) {

 if (!is.structureSim(object)) stop("Not a structureSim object")
 if (is.null(eigenSelect)) eigenSelect <- c(1:dim(object$details$eigenvalues)[2])
 
 cat("Report For a structureSim Class \n\n")
 NextMethod()
 cat(paste("Simulated eigenvalues","\n\n"))
 object$details$eigenvalues <- round(object$details$eigenvalues[,eigenSelect], ...)
 colnames(object$details$eigenvalues) <- paste("E",eigenSelect,sep="")
 print(object$details$eigenvalues)
 cat(paste("\n\n Number of factors retained by each index for each simulation","\n\n"))
 object$details$components <- round(object$details$components[,index], ...)
 print(object$details$components)
 }
 # summary(mzwick, index=1:5, eigenSelect=1:10, digits=2)
 # summary.structureSim(x)
 # summary(x)
## .................................................................

## .................................................................
print.structureSim <- function(x, index=NULL, ...) {

 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(index)) index <- c(1:dim(x$nFactors)[2])
 
 res <- x$nFactors[,index]
 print(res, ...)
 }
 # print(mzwick, index=c(1:13), 2)
 # print.structureSim(x)
 # print(x)
## .................................................................

## .................................................................
boxplot.structureSim <- function(x, nFactors=NULL, eigenSelect=NULL,
                                 vLine="green", xlab="Factors",
                                 ylab="Eigenvalues", main="Eigen Box Plot", ...) {
                                 
 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(eigenSelect)) eigenSelect <- c(1:dim(x$details$eigenvalues)[2])

 boxplot(x$details$eigenvalues[,eigenSelect], xlab=xlab, ylab=ylab, main=main, ...)
 abline(v=nFactors, lty=2, col=vLine)
 }
 # boxplot(mzwick, nFactors=3, eigenSelect=1:5, vLine="blue", col="red")
 # boxplot.structureSim(x)
 # boxplot(x)
## .................................................................

## .................................................................
plot.structureSim <- function(x, nFactors=NULL, index=NULL, main="Index Acuracy Plot", ...) {

 if (!is.structureSim(x)) stop("Not a structureSim object")
 if (is.null(index)) index <- c(1:dim(x$details$components)[2])

 if (!exists("col")  == TRUE) col  <- "black"
 ylab              <- "Average Number of Factors Retained"
 tx                <- t(x[[2]][,index])
 tx                <- data.frame(Index=rownames(tx),tx)
 colnames(tx)[2]   <- "Mean"
 tx                <- tx[order(tx[,1]),]
 plot(Mean ~ Index, type="n", data=tx, main=main, ...)
 #plot(Mean ~ Index, data=tx, cex.lab=1, cex.axis=0.7, type="n", ylab=ylab)
 abline(h=nFactors, ...)
 abline(h=median(tx[2,], na.rm=TRUE), lty=2, col="black")
 for (i in 1:length(tx[,2])) lines(y=c(0,tx[i,2]), x=c(i,i), lty=2)
 }
 # plot.structureSim(x=mzwick, nFactors=3, index=c(1:10), cex.axis=0.7, col="red")
 # plot.structureSim(x)
 # plot(x)
## .................................................................

## .................................................................
is.structureSim <- function(object) {
 if (class(object) == "structureSim") return(TRUE) else return(FALSE)
 }
 # is.structureSim(mzwick)
 # is.structureSim(x)
## .................................................................