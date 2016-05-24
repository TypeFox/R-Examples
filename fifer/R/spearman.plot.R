##' Plots the relationship between two variables using a Spearman Plot
##'
##' Often data are not normally distributed, requiring the use of a spearman correlation to determine their relationship. However, doing so
##' makes it difficult to visualize the data since scatterplots of raw data present the data as if a pearson correlation were used. This function
##' plots the ranks of the data, while plotting along the axes the distributions of the raw data.
##' @title Spearman plot
##' @param x either a matrix with two columns or a vector (if y is not \code{NULL})
##' @param y a vector
##' @param dcol the color of the lines drawn for the density plot
##' @param lhist the number of breaks in the histogram
##' @param num.dnorm the number of breaks in the density line
##' @param plot.cor logical. Should the spearman correlation be outputted in the plot?
##' @param ... arguments passed to \code{plot}
##' @export
##' @author Dustin Fife
##' @examples
##' ### generate skewed data
##' x = rnorm(1000)^2
##' y = .6*x + rnorm(1000, 0, sqrt(1-.6^2))
##'
##' spearman.plot(cbind(x,y), col="red", lhist=50)
##' spearman.plot(x=iris$Sepal.Length, y=iris$Sepal.Width)
spearman.plot <- function(x, y=NULL,dcol="blue", lhist=20, num.dnorm=5*lhist, plot.cor = TRUE,...){

	if (!is.null(y)){
		x = cbind(x,y)
	}

	if (ncol(x)!=2 & is.null(y)){
		stop("You must supply either a nx2 matrix, or a y vector")
	}
	 
    ## set up layout and graphical parameters
    layMat <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    layout(layMat, widths=c(5/7, 2/7), heights=c(2/7, 5/7))
    ospc <- 0.5 # outer space
    pext <- 2 # par extension down and to the left
    bspc <- 1 # space between scatter plot and bar plots
    par. <- par(mar=c(pext, pext, bspc, bspc),
                oma=rep(ospc, 4)) # plot parameters
 
    ## scatter plot
    par(mar=c(3, 3, 1, 1), mgp=c(1.5, .35, 0), tck=-.01, cex.axis=.8)    
    plot(apply(x, 2, rank), ...)
    abline(lm(rank(x[,2])~rank(x[,1])))    
    
   	#### make axis labels
	#axis(1, at=quantile(rank(x[,1])), labels=round(quantile(x[,1]), digits=2))
	#axis(2, at=quantile(rank(x[,2])), labels=round(quantile(x[,2]), digits=2))

    ## determine barplot and height parameter
    ## histogram (for barplot-ting the density)
    xhist <- hist(x[,1], plot=FALSE, breaks=seq(from=min(x[,1]), to=max(x[,1]),
                                     length.out=lhist))
    yhist <- hist(x[,2], plot=FALSE, breaks=seq(from=min(x[,2]), to=max(x[,2]),
                                     length.out=lhist)) # note: this uses probability=TRUE
    ## determine the plot range and all the things needed for the barplots and lines
    xx <- seq(min(x[,1]), max(x[,1]), length.out=num.dnorm) # evaluation points for the overlaid density
    xy <- dnorm(xx, mean=mean(x[,1]), sd=sd(x[,1])) # density points
    yx <- seq(min(x[,2]), max(x[,2]), length.out=num.dnorm)
    yy <- dnorm(yx, mean=mean(x[,2]), sd=sd(x[,2]))
    ## barplot and line for x (top)
    par(mar=c(1, pext, 0, 0))
    a = barplot(xhist$density, axes=FALSE, ylim=c(0, max(xhist$density, xy)),
            space=0) # barplot
    lines(seq(from=0, to=lhist-1, length.out=num.dnorm), xy, col=dcol) # line
    axis(1, at=seq(1, to=nrow(a), length.out=4)-.5, labels=round(seq(from=min(x[,1]), to=max(x[,1]), length.out=4), digits=2))
	if (plot.cor){
		corval = cor(x, use="pairwise.complete.obs", method="spearman")[1,2]
		corval = paste("r=", round(corval, digits=3), sep="")
	    legend("topright", legend=corval, bty="n")
    }
    
    ## barplot and line for y (right)
    par(mar=c(pext, 1, 0, 0), mgp=c(1, .5, 0), tck=-.01, cex.axis=.8)
    b = barplot(yhist$density, axes=FALSE, xlim=c(0, max(yhist$density, yy)),
            space=0, horiz=TRUE) # barplot
    lines(yy, seq(from=0, to=lhist-1, length.out=num.dnorm), col=dcol) # line
    axis(2, at=seq(1, to=nrow(b), length.out=4)-.5, labels=round(seq(from=min(x[,2]), to=max(x[,2]), length.out=4), digits=2))	
	## Set up x axis with tick marks alone
	# axis(1, at = seq(from=1, to=nrow(b), length.out=4), labels = FALSE, tick=T)
	# labels <- round(seq(from=min(x[,2]), to=max(x[,2]), length.out=4), digits=2)
	# text(b[1,1] - 0.51, seq(from=1, to=nrow(b), length.out=4)-1.5, 
			# srt = 90, adj = c(0,0),labels = labels, xpd = TRUE, cex=.8)


    ## restore parameters
    par(par.)
}
