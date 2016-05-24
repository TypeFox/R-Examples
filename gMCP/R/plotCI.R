#' Plot confidence intervals
#' 
#' A function for convenient plotting of confidence intervals.
#' 
#' 
#' @param ci a (named) matrix containing the lower confidence bounds in the
#' first column, the point estimates in the second and the upper confidence
#' bounds in the third column.
#' @author Code adapted from plotCII from Frank Schaarschmidt
#' @keywords hplot
#' @examples
#' 
#' 
#' est <- c("H1"=0.860382, "H2"=0.9161474, "H3"=0.9732953)
#' # Sample standard deviations:
#' ssd <- c("H1"=0.8759528, "H2"=1.291310, "H3"=0.8570892)
#' 
#' pval <- c(0.01260, 0.05154, 0.02124)/2
#' 
#' ci <- simConfint(BonferroniHolm(3), pvalues=pval, 
#'   	confint="t", df=9, estimates=est, alpha=0.025, alternative="greater")
#' 
#' plotSimCI(ci)
#' 
#' 
#' @export plotSimCI
#' 
plotSimCI <- function(ci) {	
  # Code taken (and adapted) from package MCPAN under GPL.
  # Authors: Frank Schaarschmidt, Daniel Gerhard, Martin Sill
	
	estimate = ci[,2]
	lower=ci[,1] 
	upper=ci[,3]
	
	lines=NULL
	lineslty=2
	lineslwd=1
	linescol="black"
	CIlty = 1
	CIlwd=1
	CIcex=1
	CIcol="black"
	CIlength=NULL
	
	old.par <- par(no.readonly=TRUE)
	
	aargs <- list()
	
# check input variables
	
	if(length(estimate)<1 | (!is.numeric(estimate)&!is.integer(estimate))) {
		stop("Argument estimate should be a numeric vector")
	}
	
	k<-length(estimate)
	num <- 1:k
	
	if(is.null(names(estimate))) {
		compn <- paste("C", num, sep="")
	} else {
		compn <- names(estimate)
	}
	
	if(!is.null(lower))	{
		if(!is.numeric(lower)&!is.integer(lower)) stop("Argument lower should be a numeric vector")
		if(length(lower)!=k) stop("Argument lower should be a vector of the same length as estimate!")
	}
	
	if(!is.null(upper))	{
		if(!is.numeric(upper)&!is.integer(upper))
		{stop("Argument upper should be a numeric vector")}
		if(length(upper)!=k)
		{stop("Argument upper should be a vector of the same length as estimate!")}
	}
	
	if(!is.null(lines))	{
		if(!is.numeric(lines)&!is.integer(lines))
		{stop("Argument lines should be a numeric vector")}
	}
	
	CIlty<-rep(CIlty, length.out=k)
	CIlwd<-rep(CIlwd, length.out=k)
	CIcol<-rep(CIcol, length.out=k)
	
	mymai <- par("mai")
	
# define the plot range
	
	if(is.null(lower) & is.null(upper)){warning("No confidence limits specified: arguments lower and upper are both empty!")}
	
# define the plot range
	
	allpoints <- c(lower, estimate, upper)
	if(all(!is.finite(allpoints))){stop("Arguments estimate, lower and upper contain only infinity or missing values!")}
	allexist<-allpoints[is.finite(allpoints)]
	lplot <- min(c(allexist, lines))
	uplot <- max(c(allexist, lines))
	
# Define the final plot ranges:
	
	if ((lplot-uplot)!=0) {
		llplot <- lplot - 0.1*abs(lplot-uplot)
		uuplot <- uplot + 0.1*abs(lplot-uplot)
	} else {
		llplot <- lplot -1.2
		uuplot <- uplot +1.2
	}
	
# define the type of interval drawn,
# appropriate for unbounded CI
	
	if(is.null(lower)){
		lower<-rep(llplot,k); code<-rep(2,k)
		warning("No lower limits specified!")
	} else if(is.null(upper)){
		upper<-rep(uuplot,k); code<-rep(1,k)
		warning("No upper limits specified!")
	} else {
		code<-rep(3,k)
		
		infl<-!is.finite(lower)
		lower[infl]<-llplot
		code[infl]<-2
		
		infu<-!is.finite(upper)
		upper[infu]<-uuplot
		code[infu]<-1
		
		infts <- infl&infu
		code[infts]<-0
	}
	
# Define the defaults for main, sub, ylab, xlab:
	
	aargs$main <- "" 
	aargs$sub  <- ""
	aargs$ylab <- "" 
	aargs$xlab <- ""
	
# Box arguments
	
	BTY <-"o"
	
	plot.new()
	
	# adjust margin under the x axis according to length of comparison names
	ywidth <- 1.5 * max(strwidth(compn, units = "inches", cex = par("cex.axis"))) 
	
	if (mymai[2] < ywidth) mymai[2] <- ywidth
	
	par(mai=mymai, new=TRUE)
	
	rnum<-rev(num)
	
	aargs$y <- rnum
	aargs$x <- estimate
	aargs$axes <- FALSE
	
	aargs$xlim <- c(llplot, uuplot)
	aargs$type <- "p"
	aargs$pch <- 16
	aargs$cex <- CIcex
	
	do.call("plot", aargs)
	
	axis(side = 2, at = rnum, labels=compn, las=2)
	axis(side = 1)
	box(bty=BTY)
	
	abline(h=num, col="lightgrey", lty=3)
	
	abline(v=lines, lty=lineslty, lwd=lineslwd, col=linescol)
	
	if(is.null(CIlength)) arrlength<-1/(k*2)
	
	for(i in 1:length(num))	{		
		arrows(y0=rnum[i], y1=rnum[i], x0=lower[i], x1=upper[i],
				length = arrlength, angle = 90, code = code[i],
				col = CIcol[i], lty = CIlty[i], lwd = CIlwd[i])
	}
	
	par(old.par)	
}