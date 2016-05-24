#' linear regression with plotting
#' 
#' uses \code{\link{lm}}; plots data if add=FALSE, draws the regression line
#' with \code{\link{abline}} and writes the formula with \code{\link{legend}}
#' 
#' @return None, used for plotting and drawing.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2011-2012, 2015
#' @seealso \code{\link{lm}}, \code{\link{mReg}}, \code{\link{expReg}}, \code{\link{legend}}, \code{\link{par}}, \code{\link{abline}}.
#' @keywords hplot aplot regression
#' @export
#' @examples
#' 
#' a <- 1:30
#' b <- a/2.345+rnorm(30,0,3)
#' 
#' linReg(a,b)
#' linReg(a,b, ylab="Hallo", pch=1, col=3, main="Regression by Berry")
#' linReg(a, b, pos1=15, pos2=0) # position of topleft corner of legend
#' linReg(a, b, pos1=NA, col="orange") # to suppress legend
#' 
#' # Formula specification:
#' linReg(b~a)
#' linReg(Volume~Height, data=trees)
#' 
#' # For more flexibility with the datapoints, plot first, then use linReg with add=TRUE:
#' plot(a,b, xlim=c(-5,45))
#' linReg(a, b, pos1="bottomright", add=TRUE, inset=.1) # inset: distance from plot border
#' linReg(a, b, digits=c(7,4,3), add=TRUE, col=3, lty=2, lwd=4, level=0.8)
#' linReg(a, b, pos1="topleft", inset=c(-0.1, 0.3), legargs=list(xpd=TRUE), add=TRUE)
#' 
#' @param x Numeric or formula (see examples). Vector with values of explanatory variable
#' @param y Numeric. Vector with values of dependent variable. DEFAULT: NULL
#' @param data Dataframe. If x is a formula, the according columns from data are used as x and y. DEFAULT: NULL
#' @param add Logical. If TRUE, line and text are added to the existing graphic. DEFAULT: FALSE (plots datapoints first and then the line.)
#' @param digits Numeric vector of length \eqn{\ge 1}. Specifies number of digits a,b,r,e are rounded to 
#'        in the formula "y=a*x+b \\n R^2=r \\n RMSE=e", respectively.
#'        If values are not specified, they are set equal to the first. DEFAULT: 2
#' @param pch Point Character of datapoints, see \code{\link{par}}. DEFAULT: 16
#' @param col Color of the regression line, see \code{\link{par}}. DEFAULT: 2
#' @param colband Color of the confidence region band. DEFAULT: addAlpha(col)
#' @param level Confidence level, see \code{\link{predict.lm}}. DEFAULT: 0.95
#' @param lwd Numeric. Linewidth, see \code{\link{par}}. DEFAULT: 1
#' @param xlab Axis label if add=FALSE. DEFAULT: deparse(substitute(x))
#' @param ylab Axis label if add=FALSE. DEFAULT: deparse(substitute(y))
#' @param main Title if add=FALSE. Changed (if not specified) for x=formula with data. DEFAULT: "linear regression"
#' @param pos1 \code{\link{xy.coords}}-acceptable position of the formula. DEFAULT: "top"
#' @param pos2 For numerical coordinates, this is the y-position. DEFAULT: NULL, as in \code{\link{legend}}
#' @param inset Numeric vector of length \eqn{\le 2}. inset distance(s) from the margins as a fraction of the plot region when formula legend is placed by keyword. DEFAULT: 0
#' @param legargs list of arguments passed to legend, like list(cex=0.8, xpd=TRUE, bg="white"), ...  
#'        xpd specifies whether formula can be written only inside the plot region (when FALSE) 
#'        or inside the figure region including mar (when TRUE) or in the entire device region including oma (when NA). DEFAULT: NULL
#' @param \dots Further arguments passed to \code{\link{plot}} and \code{\link{abline}}.
#' 
linReg <- function(
x,
y=NULL,
data=NULL,
add=FALSE,
digits=2,
pch=16,
col=2,
colband=addAlpha(col),
level=0.95,
lwd=1,
xlab=deparse(substitute(x)),
ylab=deparse(substitute(y)),
main="linear regression",
pos1="top",
pos2=NULL,
inset=0,
legargs=NULL,
...)
{
if(inherits(x,"formula"))
  {
  mf <- model.frame(x, data=data)
  x <- mf[,2]
  y <- mf[,1]
  if(missing(xlab)) xlab <- colnames(mf)[2]
  if(missing(ylab)) ylab <- colnames(mf)[1]
  if(!missing(data) & missing(main)) main <- paste("linear regression of",
                                                    deparse(substitute(data)))
  }
# make new plot if add is FALSE (the default):
if (!add) plot(x, y, las=1, pch=pch, xlab=xlab, ylab=ylab, main=main, ...)
# do linear regression and plotting
mod <- lm( y ~ x )
abline(mod, col=col, lwd=lwd, ...)
x2 <- seqR(par("usr")[1:2], len=100)
pred <- predict(mod, newdata=data.frame(x=x2), interval="confidence", level=level )
ciBand(yu=pred[,3], yl=pred[,2], x=x2, colb=colband, add=TRUE)
# expand digits vector, if necessary
digits <- rep(digits, length.out=4)
# Prepare formula writing
a <- round( coef(mod)[2] , digits[1])
b <- round( coef(mod)[1] , digits[2])
r <- round( summary(mod)$r.squared , digits[3])
Txt <- paste0("y = ", a, " * x ", ifelse(b>0," + ", " - "), abs(b), "\nR\U00B2 = ",
          r, "\nRMSE = ", signif(rmse(x,y), digits[4]))
# write formula
legdef <- list(x=pos1, y=pos2, legend=Txt, bty="n", text.col=col, inset=inset)
do.call(legend, owa(legdef, legargs, "legend", "inset"))
} # end of function
