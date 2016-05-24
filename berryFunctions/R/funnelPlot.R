#' Funnel plots for proportional data
#' 
#' Funnel plots for proportional data with confdence interval based on sample size. Introduced by Stephen Few, 2013
#' 
#' @return Nothing - the function just plots
#' @note the default for lty is not taken from par("lty"). This would yield "solid".
#'  Overwriting lty for one of the three line categories then produces
#'  eg c("2", "solid", "solid"), which cannot be processed by legend.\cr
#'  \bold{Wilson's Method:} algebraic approximation to the binomial distribution, very accurate, even for very small numbers.\cr
#'  \url{http://www.apho.org.uk/resource/item.aspx?RID=39445} see "contains".\cr
#'  \bold{classic = Stephen Few's Method = the way I knew it:} sqrt( mu*(1-mu) / n )\cr
#'  \url{http://www.jerrydallal.com/LHSP/psd.htm}\cr
#'  \url{http://commons.wikimedia.org/wiki/File:ComparisonConfidenceIntervals.png}\cr
#'  The apho Wilson method first yielded wrong upper limits in my translation (it needs 0:1 instead of \%). Thus I added the wikipedia formula:\cr
#'  \url{http://de.wikipedia.org/wiki/Konfidenzintervall_einer_unbekannten_Wahrscheinlichkeit#Wilson-Intervall}\cr
#'  \url{http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval}\cr
#'  Which other methods should I include? (That's not the hard part anymore)
#' @section The basic idea: Salesman A (new to the job) has had 3 customers and
#' sold 1 car. So his success rate is 0.33. Salesman B sold 1372 customers 632
#' cars, thus having a success rate of 0.46 Promoting B solely because of the
#' higher rate fails to take experience and opportunity (n) into account! This
#' dilemma is what the funnel plot with the confidence interval (ci) solves.
#' See Stephen Few and Katherine Rowel's PDF for details on the interpretation.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Oct 2013
#' @references
#'    http://www.perceptualedge.com/articles/visual_business_intelligence/variation_and_its_discontents.pdf\cr
#'    \url{http://sfew.websitetoolbox.com/post/variation-and-its-discontents-6555336?}\cr
#'    Excellent explanation of bayesian take on proportions: \url{http://varianceexplained.org/r/empirical_bayes_baseball/}
#' @keywords hplot aplot
#' @export
#' @examples
#' 
#' # Taken directly from Stephen Few's PDF:
#' funnel <- read.table(header=TRUE, text="
#' Name SampleSize Incidents
#' Tony 2 2
#' Mike 400 224
#' Jan 100 54
#' Bob 1000 505
#' Sheila 2 1
#' Jeff 10 5
#' Sandy 500 236
#' Mitch 200 92
#' Mary 10 3
#' John 2 0")
#' 
#' str(funnel)
#' X <- funnel$Incidents
#' N <- funnel$SampleSize
#' 
#' barplot(X/N, names=funnel$Name, main="success rate")
#' # not showing n!
#' 
#' funnelPlot(X,N)
#' # arguments for subfunctions as text may be given this way:
#' funnelPlot(x=X, n=N, labels=funnel$Name, at=list(cex=0.7, col="red"))
#' # Labeling many points is not very clear...
#' 
#' # Even though Jan is more successfull than Mary in succes rate terms, both are
#' # easily within random variation. Mary may just have had a bad start.
#' # That Mike is doing better than average is not random, but (with 95% confidence)
#' # actually due to him being a very good seller.
#' 
#' # one more interesting option:
#' funnelPlot(X,N, a3=list(lty=2))
#' 
#' funnelPlot(X,N, a3=list(col=2, lwd=5))
#' # changing round line ends in legend _and_ plot is easiest with
#' par(lend=1)
#' funnelPlot(X,N, a3=list(col=2, lwd=5))
#' 
#' # The Wilson method yields slighty different (supposedly better) limits for small n:
#' funnelPlot(X,N, method="classic", al=list(title="Standard Method"))
#' funnelPlot(X,N, add=TRUE, method="wilson", a3=list(lty=2, col="red"),
#'            a2=list(lty=2, col="blue"), al=list(x="bottomright", title="Wilson Method"))
#' 
#' # Both Wilson method implementations yield the same result:
#' funnelPlot(X,N, method="wilson")
#' funnelPlot(X,N, add=TRUE, method="wilsonapho",
#'            a3=list(lty=2, col="red"), a2=list(lty=2, col="blue"))
#' 
#' 
#' # Note on nl used in the function, the n values for the ci lines:
#' plot(     seq(      10 ,       300 , len=50), rep(  1, 50) )
#' points(10^seq(log10(10), log10(300), len=50), rep(0.8, 50) )
#' abline(v=10)
#' # CI values change rapidly at small n, then later slowly.
#' # more x-resolution is needed in the first region, so it gets more of the points
#' 
#' @param x Numeric vector with number of successes (cases).
#' @param n Numeric vector with number of trials (population).
#' @param labels Labels for points. DEFAULT: NULL
#' @param method Method to calculate Confidence interval, see "note" below. Can also be "wilson". DEFAULT: "classic"
#' @param add Add to existing plot instead of drawing new plot? DEFAULT: FALSE
#' @param xlim Graphical parameters, see \code{\link{par}} and \code{\link{plot}}. DEFAULT: range(n, finite=TRUE)
#' @param ylim y limit in [0:1] DEFAULT: range(x/n*100, finite=TRUE)
#' @param las DEFAULT: 1
#' @param xlab DEFAULT: "Sample size n"
#' @param ylab DEFAULT: "Success rate [\%]"
#' @param main DEFAULT: "Funnel plot for Proportions"
#' @param a3 List with arguments for CI lines at 3*sd (eg: col, lty, lwd, lend, etc.). 
#'        Overwrites defaults that are defined within the function (if contentually possible). DEFAULT: NULL
#' @param a2 Arguments for line of 2 sd. DEFAULT: NULL
#' @param am Arguments for mean line. DEFAULT: NULL
#' @param ap Arguments for the data points (cex, etc.). DEFAULT: NULL
#' @param at Arguments for text (labels of each point). DEFAULT: NULL
#' @param al Arguments for \code{\link{legend}} (text.col, bty, border, y.intersp, etc.). DEFAULT: NULL
#' @param \dots further arguments passed to plot only!
#' 
funnelPlot <- function(
x,
n,
labels=NULL,
method="classic",
add=FALSE,
xlim=range(n, finite=TRUE),
ylim=range(x/n*100, finite=TRUE),
las=1,
xlab="Sample size n",
ylab="Success rate [%]",
main="Funnel plot for Proportions",
a3=NULL,
a2=NULL,
am=NULL,
ap=NULL,
at=NULL,
al=NULL,
...)
{
# Data (proportions) -----------------------------------------------------------
p <- x/n*100 # p: proportion of success
m <- mean(p, na.rm=TRUE) # m: mean value of proportions
# Distribute line point values with high density at rapidly changing curves
nl <- 10^seq(log10(xlim[1]), log10(xlim[2]), len=500) # nl: n for ci-lines
# calculate CI  (2*sd and 3*sd limits) -----------------------------------------
# f: factor for confidence interval
f1 <- qnorm(1-0.05/2) # 1.959964   1.96    # 0.025 = 2.5%    alpha=0.95
f2 <- qnorm(1-0.002/2) # 3.090232          # 0.001 = 0.1%    alpha=0.998
if(method=="wilsonapho")
ci <- data.frame(l2sigma=(2*nl*m/100+f1^2 - f1*sqrt(f1^2+4*nl*m/100*(1-m/100)))/(nl+f1^2)/2,
                 u2sigma=(2*nl*m/100+f1^2 + f1*sqrt(f1^2+4*nl*m/100*(1-m/100)))/(nl+f1^2)/2,
                 l3sigma=(2*nl*m/100+f2^2 - f2*sqrt(f2^2+4*nl*m/100*(1-m/100)))/(nl+f2^2)/2,
                 u3sigma=(2*nl*m/100+f2^2 + f2*sqrt(f2^2+4*nl*m/100*(1-m/100)))/(nl+f2^2)/2 )*100
else
if(method=="wilson")
ci <- data.frame(l2sigma=1/(1+f1^2/nl)*(m/100+f1^2/2/nl - f1*sqrt(m/100*(1-m/100)/nl+f1^2/4/nl^2)),
                 u2sigma=1/(1+f1^2/nl)*(m/100+f1^2/2/nl + f1*sqrt(m/100*(1-m/100)/nl+f1^2/4/nl^2)),
                 l3sigma=1/(1+f2^2/nl)*(m/100+f2^2/2/nl - f2*sqrt(m/100*(1-m/100)/nl+f2^2/4/nl^2)),
                 u3sigma=1/(1+f2^2/nl)*(m/100+f2^2/2/nl + f2*sqrt(m/100*(1-m/100)/nl+f2^2/4/nl^2)) )*100
else
if(method=="classic")
ci <- data.frame(l2sigma= m + f1*sqrt( m*(100-m) / nl ),
                 u2sigma= m - f1*sqrt( m*(100-m) / nl ),
                 l3sigma= m + f2*sqrt( m*(100-m) / nl ),
                 u3sigma= m - f2*sqrt( m*(100-m) / nl ))
else stop("Wrong method. Possible are 'classic' and 'wilson'")
#
# Plot preparation -------------------------------------------------------------
# default arguments (da)
da3 <- list(x=nl, col="black", lty=1, lwd=par("lwd")) # da_: default arguments for _lines at 3*sd
da2 <- list(x=nl, col="gray",  lty=1, lwd=par("lwd")) # lines 2sd
dam <- list(h=m, col="orange", lty=1, lwd=par("lwd")) # abline mean
dap <- list(x=n, y=p, col="orange", pch=16)  # points
dat <- list(x=n, y=p, labels=labels, adj=c(0,1), col="black") # text (labels)
# final arguments for lines (so they can also be used by legend default arguments)
fa3 <- owa(da3, a3, "y")
fa2 <- owa(da2, a2, "y")
fam <- owa(dam, am, "h")
# default arguments legend:
dal <- list(x="topright", legend=c("3*sd CI (99.8%)","2*sd CI (95%)", "mean"),
            col=c(fa3$col, fa2$col, fam$col),
            lty=c(fa3$lty, fa2$lty, fam$lty),
            lwd=c(fa3$lwd, fa2$lwd, fam$lwd))
#browser()
# Plot everything --------------------------------------------------------------
if(!add) plot(nl, ci[,1], type="n", ylim=ylim, las=las, xlab=xlab, ylab=ylab, main=main, ...)
do.call(abline, args = fam ) # fam: final arguments for mean line drawing
do.call(lines,  args = c(list(y=ci[,1]), fa2) )
do.call(lines,  args = c(list(y=ci[,2]), fa2) )
do.call(lines,  args = c(list(y=ci[,3]), fa3) )
do.call(lines,  args = c(list(y=ci[,4]), fa3) )
do.call(points, args = owa(dap, ap, "x","y") )
do.call(text,   args = owa(dat, at) )
do.call(legend, args = owa(dal, al, "col", "lty", "lwd", "pch") )
} # End of function ------------------------------------------------------------
