#' Exponential regression with plotting
#' 
#' uses \code{\link{lm}}; plots data if add=FALSE, draws the regression line
#' with \code{\link{abline}} and confidence interval with \code{\link{polygon}}
#' and writes the formula with \code{\link{legend}}
#' 
#' @return \code{\link{predict.lm}} result.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec. 2014
#' @seealso \code{\link{lm}}, \code{\link{mReg}}, \code{\link{linReg}}.
#' @keywords hplot aplot regression
#' @export
#' @examples
#' 
#' x <- runif(100, 1, 10)
#' y <- 10^(0.3*x+rnorm(100, sd=0.3)+4)
#' plot(x,y)
#' expReg(x,y)
#' expReg(x,y, logy=FALSE)
#' expReg(x,y, predictnew=6, plot=FALSE)
#' expReg(x,y, predictnew=3:6, interval="none", plot=FALSE)
#' 
#' @param x Numeric or formula (see examples). Vector with values of explanatory variable
#' @param y Numeric. Vector with values of dependent variable. DEFAULT: NULL
#' @param data Dataframe. If x is a formula, the according columns from data are used as x and y. DEFAULT: NULL
#' @param logy Plot with a logarithmic y axis?  Calls \code{\link{logAxis}}. DEFAULT: TRUE
#' @param predictnew Vector with values to predict outcome for. Passed as \code{newdata} to \code{\link{predict.lm}}. DEFAULT: NULL
#' @param interval Interval for prediction. DEFAULT: "confidence"
#' @param plot Plot things at all? If FALSE, predictnew will still be returned. DEFAULT: TRUE
#' @param digits Numeric vector of length \eqn{\ge 1}. Specifies number of digits a,b,r,e are rounded to 
#'        in the formula "y=a*log(x)+b, R^2, RMSE=e", respectively. 
#'        If values are not specified, they are set equal to the first. DEFAULT: 2
#' @param inset Numeric vector of length \eqn{\le 2}. inset distance(s) from the margins 
#'        as a fraction of the plot region when formula is placed by keyword. DEFAULT: 0
#' @param xpd Logical, specifying wheter formula can be written only inside the plot region (when FALSE) 
#'        or inside the figure region including mar (when TRUE)
#'        or in the entire device region including oma (when NA). DEFAULT: par("xpd")
#' @param pos1 \code{\link{xy.coords}}-acceptable position of the formula. DEFAULT: "top"
#' @param pos2 For numerical coordinates, this is the y-position. DEFAULT: NULL, as in \code{\link{legend}}
#' @param add Logical. If TRUE, line and text are added to the existing graphic. DEFAULT: FALSE (plots datapoints first and then the line.)
#' @param pch Point Character, see \code{\link{par}}. DEFAULT: 16
#' @param col Color of points, see \code{\link{par}}. DEFAULT: rgb(0,0,0, 0.5)
#' @param modcol color of model line. DEFAULT: 2
#' @param lwd Numeric. Linewidth, see \code{\link{par}}. DEFAULT: 1
#' @param xlab,ylab,main Character / Expression. axis label and graph title if add=FAlSE. DEFAULT: internal from names
#' @param xlim,ylim graphic range. DEFAULT: range(x)
#' @param \dots Further arguments passed to \code{\link{plot}} and \code{\link{abline}}.
#' 
expReg <- function(
x,
y=NULL,
data=NULL,
logy=TRUE,
predictnew=NULL,
interval="confidence",
plot=TRUE,
digits=2,
inset=0,
xpd=par("xpd"),
pos1="top",
pos2=NULL,
add=FALSE,
pch=16,
col=rgb(0,0,0, 0.5),
modcol=2,
lwd=1,
xlab=deparse(substitute(x)),
ylab=deparse(substitute(y)),
main="exponential regression",
xlim=range(x),
ylim=range(y),
...)
{
# deparse labs before x and y are evaluated:
xlab <- xlab ;  ylab <- ylab
# Formula:
if(inherits(x,"formula"))
  {
  if(!missing(data))
    {                   # get x and y from data.frame
    name <- as.character(x)[-1]
    x <- data[ , name[2] ]  ;  if(missing(xlab)) xlab <- name[2]
    y <- data[ , name[1] ]  ;  if(missing(ylab)) ylab <- name[1]
    if(missing(main)) main <- paste("exponential regression of",
                                    deparse(substitute(data)))
    } else
    {                   # get x and y from formula directly
    name <- as.character(x)[-1]
    x <- get(name[2]) ;  if(missing(xlab)) xlab <- name[2]
    y <- get(name[1]) ;  if(missing(ylab)) ylab <- name[1]
    }
  }
if(is.null(y)) stop("y cannot be NULL. formula statement might be wrong.")
x <- x[y!=0]
y <- y[y!=0]
y <- log10(y)
# Model:
mod <- lm( y ~ x )
# expand digits vector, if necessary
digits <- rep(digits, length.out=4)
# Prepare formula writing
a <- round( coef(mod)[2] , digits[1])
b <- round( coef(mod)[1] , digits[2])
Txt <- paste0("y = 10^(", a, " * x ", ifelse(b>0," + ", " - "), abs(b), ")\nR\U00B2 = ",
          signif(rsquare(x,y), digits[3]), "\nRMSE = ", signif(rmse(x,y), digits[4]))
# Acutal plotting:
if(!logy) y <- 10^y
if(plot){
# make new plot if add is FALSE (the default):
if (!add) plot(x, y, pch=pch, xlab=xlab, ylab=ylab, main=main, yaxt="n",
               xlim=xlim, ylim=ylim, if(logy) type="n", col=col, ...)
# logAxis:
if(logy)
  {
  logAxis(2)
  points(x, y, las=1, pch=pch, col=col, ...)
  } else axis(2, las=1)
# Model line
lx <- seq(par("usr")[1], par("usr")[2], length.out=100)
ly <- predict(mod, newdata=data.frame(x=lx), interval="confidence" )
if(!logy) ly <- 10^ly
# confidence interval band
polygon(x=c(lx, rev(lx)), y=c(ly[,2], rev(ly[,3])), col=rgb(0.3,0.5,0, 0.5))
# prediction line
lines(lx, ly[,1], col=modcol, lwd=lwd, ...)
# write formula
legend(pos1, pos2, Txt, bty="n", text.col=modcol, inset=inset, xpd=xpd)
} # end if plot
# return output if wanted:
if(!is.null(predictnew))
   10^as.numeric(predict(mod, newdata=data.frame(x=predictnew), interval=interval))
} # end of function
