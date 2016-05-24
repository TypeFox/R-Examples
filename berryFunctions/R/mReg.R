# multiple Regression with several function types
#  1 linear                 a*x + b
#  2 quadratic (parabola)   a*x^2 + b*x + c
#  3 kubic                  a*x^3 + b*x^2 + c*x + d
#  4 Polynom 4th degree     a*x^4 + b*x^3 + c*x^2 + d*x + e
#  5 Polynom 5              a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f
#  6 logarithmic            a*log(x) + b
#  7 exponential            a*e^(b*x)
#  8 power/root             a*x^b
#  9 reciprocal             a/x + b
# 10 rational               1 / (a*x + b)
# 11 exponential 4 Param    a*e^(b*(x+c)) + d
# Yet to add: 12 hyperbolic             sinh(), cosh(), tanh()
#             13 powerplus a^x + b

#' Multiple regression
#' 
#' Multiple regression fitting various function types including e.g. linear, cubic, logarithmic, exponential, power, reciprocal. 
#' Quick way to find out what function type fits the data best. 
#' Plots data and fitted functions and adds a legend with the functions (or their types=structure) sorted by R squared. 
#' Returns the fitted functions with their parameters and R^2 values in a data.frame.
#' 
#' @details legendform : example\cr 
#'          full : 7.8*x + 6.31\cr 
#'          form : a*x+b\cr 
#'          nameform : linear a*x+b\cr 
#'          name : linear\cr\cr 
#'          full can be quite long, especially with Poly45=TRUE!
#' 
#' @return data.frame with rounded R squared, formulas, and full R^2 and parameters for further use. 
#'         Rownames are the names (types) of function. Sorted decreasingly by R^2
#' @note If you're adjusting the appearance (lwd, lty, col) of single lines,
#' set parameters in the following order:\cr 
#' # 1 linear a*x + b\cr 
#' # 2 quadratic (parabola) a*x^2 + b*x + c\cr 
#' # 3 kubic a*x^3 + b*x^2 + c*x + d\cr
#' # 4 Polynom 4th degree a*x^4 + b*x^3 + c*x^2 + d*x + e\cr 
#' # 5 Polynom 5 a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f \cr 
#' # 6 logarithmic a*log(x) + b \cr
#' # 7 exponential a*e^(b*x) \cr 
#' # 8 power/root a*x^b \cr 
#' # 9 reciprocal a/x + b \cr 
#' # 10 rational 1 / (a*x + b) \cr 
#' # 11 exponential 4 Param a*e^(b*(x+c)) + d \cr
#' 
#' Negative values are not used for regressions containing logarithms; with warning.\cr
#' exp_4par was originally developed for exponential temperature decline in a cup of hot water.
#' @section warning: A well fitting function does NOT imply correct causation!\cr 
#'         A good fit does NOT mean that you describe the bahaviour of a system adequatly!\cr 
#'         Extrapolation can be DANGEROUS!\cr 
#'         Always extrapolate to see if a function fits the expected results there as well.\cr 
#'         Avoid overfitting: Poly45 will often yield good results (in terms of R^2), but can be way overfitted. 
#'         And outside the range of values, they act wildly.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, Dec 2012, updated April and Aug 2013, sept 2015
#' @seealso \code{\link{glm}}, \code{\link{lm}}, \code{\link{optim}}
#' @references Listed here: \url{http://rclickhandbuch.wordpress.com/berryfunctions/#mReg}
#' @keywords aplot hplot regression nonlinear multivariate
#' @export
#' @examples
#' 
#' set.seed(12)
#' x <- c(runif(100,0,3), runif(200, 3, 25)) # random from uniform distribution
#' y <- 12.367*log10(x)+7.603+rnorm(300)     # random from normal distribution
#' plot(x,y, xlim=c(0,40))
#' mReg(x,y) # warning comes from negative y-values (suppress with quiet=TRUE)
#' 
#' # Formula specification:
#' mReg(Volume~Height, data=trees)
#' 
#' # NA management
#' x[3:20] <- NA
#' mReg(x,y)
#' 
#' # Passing arguments to legend:
#' mReg(x,y, pch=1, legargs=list(x="bottomright", cex=0.7), legendform="form")
#' 
#' mReg(x,y, col=rainbow2(11))
#' mReg(x,y, extend=0.2) # less empty space around data points
#' mReg(x,y, nbest=4) # only 4 distributions plotted
#' mReg(x,y, legargs=list(x=7, y=8, bty="o", cex=0.6)) # Legend position as coordinates
#' 
#' \dontrun{
#' ## Rcmd check --as-cran doesn't like to open external devices,
#' ## so this View example is excluded from running in the checks.
#' View(mReg(x,y, Poly45=TRUE, exp_4=TRUE, plot=FALSE)) # exp_4: fit more distributions
#' }
#' # optim methods often yield different results, so be careful using this.
#' # I might insert a possibility to specify initial values for optim.
#' # 4 Parameters allow several combinations to yield similarly good results!
#' plot( 0:10, 3.5*exp(0.8*( 0:10 + 2      )) + 15 , type="l")
#' lines(0:10,  18*exp(0.8*( 0:10 - 2.5e-05)) - 5, col=2)
#' 
#' 
#' # okay, different dataset:
#' x <- c(1.3, 1.6, 2.1, 2.9, 4.4, 5.7, 6.6, 8.3, 8.6, 9.5)
#' y <- c(8.6, 7.9, 6.6, 5.6, 4.3, 3.7, 3.2, 2.5, 2.5, 2.2)
#' mReg(x,y, legargs=list(cex=0.7, x="topright"), main="dangers of extrapolation")
#' points(x,y, cex=2, lwd=2)
#' # Polynomial fits are good within the data range, but, in this case obviously,
#' # be really careful extrapolating! If you know that further data will also be low,
#' # add another point to test differences:
#' mReg(c(x,11,13,15), c(y,2,2,2), xf="myX", yf="myY", Poly45=TRUE, legendform="name")
#' points(x,y, cex=2, lwd=2)
#' # The Polynomials are still very good: they have 5 to 6 Parameters, after all!
#' # Poly45 is set to FALSE by default to avoid such overfitting.
#'  
#' mReg(x,y, pcol=8, ncol=0) # no return to console
#' 
#' # only plot a subset: best n fits, minimum fit quality, or user selection
#' mReg(x,y, pcol=8, ncol=2, nbest=4)
#' mReg(x,y, pcol=8, ncol=2, R2min=0.7)
#' mReg(x,y, pcol=8, ncol=2, selection=c(2,5,8))
#' # selecting the fifth degree polynomial activates Poly45 (in the output table)
#' 
#' # Add to axisting plot:
#' plot(x,y, xlim=c(0,40))
#' mReg(x,y, add=TRUE, lwd=12:1/2, ncol=0)
#' # lwd, lty can be vectors of length 12, specifying each line separately.
#' # Give those in fix order (see section notes), not in best-fit order of the legend.
#' # The order is Polynomial(1:5), log, exp, power, reciprocal, rational, exp_4_param
#' # color has to be a vector of 12
#' # opposedly, lwd and lty are repeated 12 times, if only one value is given
#' 
#' 
#' # One more dataset:
#' j <- c(5,8,10,9,13,6,2) ; k <- c(567,543,587,601,596,533,512)
#' # Inset from margin of plot region:
#' mReg(j,k, legargs=list(x="bottomright", inset=.05, bty="o"), legendform="name")
#' # Legend forms
#' mReg(j,k, legargs=list(x="bottomright"), legendform="name")
#' mReg(j,k, legargs=list(x="bottomright"), legendform="form")
#' mReg(j,k, legargs=list(x="bottomright"), legendform="nameform")
#' mReg(j,k, legargs=list(x="bottomright"), legendform="full")
#' 
#' 
#' # The question that got me started on this whole function...
#' # exponential decline of temperature of a mug of hot chocolate
#' tfile <- system.file("extdata/Temp.txt", package="berryFunctions")
#' temp <- read.table(tfile, header=TRUE, dec=",")
#' head(temp)
#' plot(temp)
#' temp <- temp[-20,] # missing value - rmse would complain about it
#' 
#' x <- temp$Minuten
#' y <- temp$Temp 
#' mReg(x,y, exp_4=TRUE, selection=11)
#' # y=49*e^(-0.031*(x - 0  )) + 25 correct, judged from the model:
#' # Temp=T0 - Te *exp(k*t) + Te     with    T0=73.76,  Tend=26.21, k=-0.031
#' # optmethod="Nelder-Mead"  # y=52*e^(-0.031*(x + 3.4)) + 26 wrong
#' 
#' 
#' x <- seq(1, 1000, 1)
#' y <- (x+22)/(x+123) # can't find an analytical solution so far. Want to check out nls
#' mReg(x, y, legargs=list(x="right"))
#' 
#' 
#' # Solitaire Results. According to en.wikipedia.org/wiki/Klondike_(solitaire):
#' # Points=700000/Time + Score
#' # I recorded my results as an excuse to play this game a lot.
#' sfile <- system.file("extdata/solitaire.txt", package="berryFunctions")
#' solitaire <- read.table(sfile, header=TRUE)
#' mReg(solitaire$Time, solitaire$Points) # and yes, reciprocal ranks highest! Play Fast!
#' mReg(solitaire$Time, solitaire$Bonus, xlim=c(50,200), extend=0, nbest=3)
#' sol <- unique(na.omit(solitaire[c("Time","Bonus")]))
#' sol
#' sol$official <- round(700000/sol$Time/5)*5
#' mReg(sol$Time, sol$Bonus, extend=0, selection=9, col=rep(4,10), legendform="full")
#' plot(sol$Time, sol$official-sol$Bonus, type="l")
#' 
#' # multivariate regression should be added, too:
#' sfile <- system.file("extdata/gelman_equation_search.txt", package="berryFunctions")
#' mv <- read.table(sfile, header=TRUE)
#' 
#' sfile <- system.file("extdata/mRegProblem.txt", package="berryFunctions")
#' x <- read.table(sfile, header=TRUE)$x
#' y <- read.table(sfile, header=TRUE)$y
#' mReg(x,y,  digits=6) # all very equal
#' x2 <- x-min(x)
#' mReg(x2,y, digits=6)          #  Formulas are wrong if digits is too low!! 
#' #mReg(x2,y, legendform="full")
#'
#' # Zero and NA testing (to be moved to unit testing someday...) 
#' mReg(1:10, rep(0,10))
#' mReg(1:10, c(rep(0,9),NA))
#' mReg(1:10, rep(NA,10))
#' mReg(rep(1,10), 1:10)
#' mReg(rep(0,10), 1:10)
#' mReg(c(rep(0,9),NA), 1:10)
#' mReg(rep(NA,10), 1:10)
#'
#'mReg(1:10, rep(0,10), quiet=TRUE)
#'mReg(1:10, c(rep(0,9),NA), quiet=TRUE)
#'mReg(1:10, rep(NA,10), quiet=TRUE)
#'mReg(rep(1,10), 1:10, quiet=TRUE)
#'mReg(rep(0,10), 1:10, quiet=TRUE)
#'mReg(c(rep(0,9),NA), 1:10, quiet=TRUE)
#'mReg(rep(NA,10), 1:10, quiet=TRUE)
#' 
#' @param x Vector with x coordinates or formula (like y~x), the latter is passed to \code{\link{model.frame}}
#' @param y Vector with y values. DEFAULT: NULL (to enable x to be a formula)
#' @param data data.frame in which formula is applied. DEFAULT: NULL
#' @param Poly45 Logical. Should 4th and 5th degree polynomials also be fitted? DEFAULT: FALSE, as the formulas are very long.
#' @param exp_4 Logical. Return 4-parametric exponential distibution fits (via \code{\link{exp4p}}) in the output table? (only best fit is plotted).
#'        exp_4par_ini has the initial values of exponential fitting with the data relocated to first quadrant. 
#'        The others are optimized with the methods of \code{\link{optim}}. DEFAULT: FALSE
#' @param xf Character. x name for Formula. DEFAULT: substitute(x) before replacing zeros in x and y
#' @param yf Ditto for y
#' @param ncolumns Number of columns in output. Set lower to avoid overcrowding the console. DEFAULT: 9
#' @param plot Logical. plot data and fitted functions? DEFAULT: TRUE
#' @param add Logical. add lines to existing plot? DEFAULT: FALSE
#' @param nbest Integer. Number of best fitting functions to be plotted (console output table always has all). DEFAULT: 12
#' @param R2min Numerical. Minimum Rsquared value for function type to be plotted. 
#'        Suggestion: 0.6 (2/3 of variation of y is explained by function of x). DEFAULT: empty
#' @param selection Integers of functions to be plotted, assigned as in list in section "note". DEFAULT: NULL, meaning all
#' @param digits Integer. number of significant digits used for rounding formula parameters and R^2 displayed. DEFAULT: 2
#' @param extend Numerical. Extention of axis ranges (proportion of range). DEFAULT: 0.4
#' @param xlim Numerical vector with two values, defining the x-range of the lines to be plotted.  DEFAULT: extended range(x)
#' @param ylim Ditto for Y-axis
#' @param xlab Character. default labels for axis labeling and for formulas. DEFAULT: substitute(x) before replacing zeros in x and y
#' @param ylab Ditto for y axis.
#' @param las Integer in 0:4. label axis style. See \code{\link{par}}. DEFAULT: 1
#' @param lwd Numerical of length 12. line width for lines. DEFAULT: rep(1,12)
#' @param lty Numerical of length 12. line type. DEFAULT: rep(1,12)
#' @param col Numerical of length 12. line colors. DEFAULT: NULL, means they are specified internally
#' @param pcol Color used for the data-points themselves. DEFAULT: par('col')
#' @param pch Integer or single character. Point CHaracter for the data points. See \code{\link{par}}. DEFAULT: 16
#' @param legend Logical. Add legend to plot? DEFAULT: TRUE
#' @param legargs List. List of arguments passed to \code{\link{legend}}. Will overwrite internal defaults. DEFAULT: NULL
#' @param legendform One of 'full', 'form', 'nameform' or 'name'. Complexity (and length) of legend in plot. See Details. DEFAULT: 'nameform'
#' @param quiet Suppress warnings about value removal (NAs, smaller 0, etc)? DEFAULT: FALSE
#' @param \dots Further graphical parameters passed to plot
#' 
mReg <- function(
x,
y=NULL,
data=NULL,
Poly45=FALSE,
exp_4=FALSE,
xf=deparse(substitute(x)),
yf=deparse(substitute(y)),
ncolumns=9,
plot=TRUE,
add=FALSE,
nbest=12,
R2min,
selection=NULL,
digits=2,
extend=0.4,
xlim=extendrange(x, f=extend), # range(x, finite=TRUE) + c(-1,1)*extend*diff(range(x, finite=TRUE))
ylim=extendrange(y, f=extend),
xlab=xf,
ylab=yf,
las=1,
lwd=rep(1,12),
lty=rep(1,12),
col=NULL,
pcol=par("col"),
pch=16,
legend=TRUE,
legargs=NULL,
legendform="nameform",
quiet=FALSE,
...)
{
# Function start
# input checking
if( (xf %in% letters[1:6] | yf %in% letters[1:6])  &  legendform %in% c("nameform", "form")  )
   if(!quiet) warning("Using single letters a to f for input variable names is not recommended, as formula forms will be difficult to read" )
if( xf=="e" | yf=="e" )
   if(!quiet) warning("Using 'e' for input variable name is not recommended, as exponential formula forms will be difficult to read" )
if(any(4:5 %in% selection)) Poly45 <- TRUE
if(11 %in% selection) exp_4 <- TRUE
if( ! round(nbest,1) %in% 0:12) stop("nbest has to be an integer between 0 and 12")
lwd <- rep(lwd, length=12)
lty <- rep(lty, length=12)
if(inherits(x,"formula"))
{
  mf <- model.frame(x, data=data)
  x <- mf[,2]
  y <- mf[,1]
  if(missing(xlab)) xf <- colnames(mf)[2]
  if(missing(ylab)) yf <- colnames(mf)[1]
  #if(!missing(data) & missing(main)) main <- paste("multiple regression of",deparse(substitute(data)))
}
# NA removal
if(any(is.na(x)|is.na(y)))
  {
  Na <- which(is.na(x)|is.na(y))
  if(!quiet) warning(length(Na), " NAs were omitted from ", length(x), " data points (",
          round(length(Na)/length(x)*100,1),"%).")
  x <- x[-Na] ; y <- y[-Na]
  } # end if NA
# vector length check:
if(length(x)!=length(y)) stop("x (",length(x), " elements) and y (",length(y),") must be of the same length.")
# log regression vectors without zeros and negative values:
neg <- which(x<=0|y<=0)
if(length(neg)!=0) 
  {
  if(!quiet) warning("For log/exp/power regressions, ",length(neg),
    " nonpositive values were removed from x and y (",round(length(neg)/length(x)*100,1),"%).")
  xg0 <- x[-neg] # xg0: x greater zero
  yg0 <- y[-neg]
  } else
  {xg0 <- x; yg0 <- y}
#
# Functions needed for function descriptions
# abbreviate parameters of fitted functions:
ab1 <- function(input) signif(input,digits)
# abbreviate parameters of fitted functions with algebraic sign (Vorzeichen):
ab <- function(input) paste0(ifelse(input>0, " + ", " - "),
                             abs(signif(input,digits)))
# Prepare Output Table
output <- as.data.frame(matrix(NA, ncol=if(Poly45) 10 else 8, nrow=11 ))
colnames(output) <- c("nr","R2","Formulas","R2full", letters[1:(ncol(output)-4)] )
#
#  1 linear --------------- a*x + b --------------------------------------------
if(length(x)>0 & length(y)>0)
  {
mod1 <- lm( y ~ x )
output$R2[1] <- summary(mod1)$r.squared
output$a [1] <- coef(mod1)[2]
output$b [1] <- coef(mod1)[1]
#  2 quadratic (parabola) - a*x^2 + b*x + c ------------------------------------
mod2 <- lm(y ~ I(x^2) + x)
output$a [2] <- coef(mod2)[2]
output$b [2] <- coef(mod2)[3]
output$c [2] <- coef(mod2)[1]
output$R2[2] <- rsquare(y, output$a[2]*x^2 + output$b[2]*x + output$c[2], quiet=quiet)
#  3 cubic ---------------- a*x^3 + b*x^2 + c*x + d ----------------------------
mod3 <- lm(y ~  poly(x,3, raw=TRUE))
output[3,5:8] <- rev(coef(mod3))
output$R2[3] <- rsquare(y, output$a[3]*x^3 + output$b[3]*x^2 + output$c[3]*x + output$d[3], quiet=quiet)
if(Poly45){
  #  4 Polynom4 ----------- a*x^4 + b*x^3 + c*x^2 + d*x + e --------------------
  mod4 <- lm(y ~  poly(x,4, raw=TRUE))
  output[4, 5:9] <- rev(coef(mod4))
  output$R2[4] <- rsquare(y, output$a[4]*x^4 + output$b[4]*x^3 + output$c[4]*x^2 + output$d[4]*x + output$e[4], quiet=quiet)
  #  5 Polynom5 ----------- a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f ------------
  mod5 <- lm(y ~  poly(x,5, raw=TRUE))
  output[5,5:10] <- rev(coef(mod5))
  output$R2[5] <- rsquare(y, output$a[5]*x^5 + output$b[5]*x^4 + output$c[5]*x^3 + output$d[5]*x^2 + output$e[5]*x + output$f[5], quiet=quiet)
  } # if Poly45 end
  } else # length(y) end
  if(!quiet) warning("linear/square/cubic not fitted (no non-NA values in dataset).")
#  6 logarithmic ---------- a*log(x) + b ---------------------------------------
if(length(xg0)>0 & length(yg0)>0)
  {
  mod6 <- lm( yg0 ~ log10(xg0) )
  output$a [6] <- coef(mod6)[2]
  output$b [6] <- coef(mod6)[1]
  output$R2[6] <- rsquare(yg0, output$a[6]*log10(xg0) + output$b[6], quiet=quiet)
#  7 exponential ---------- a*e^(b*x) ------------------------------------------
  mod7 <- lm( log(yg0) ~ xg0 )               # y = a*e^(b*x)
  output$a [7] <- exp(coef(mod7)[1])         # ln(y) = ln(a) + ln( e^(b*x) )
  output$b [7] <- coef(mod7)[2]              # ln(y) = ln(a) + b*x
  output$R2[7] <- rsquare(yg0, output$a[7]*exp(output$b[7]*xg0), quiet=quiet)
#  8 power/root ----------- a*x^b ----------------------------------------------
  mod8 <- lm( log(yg0) ~ log(xg0) )                  # y = a*x^b
  output$a [8] <- exp(coef(mod8)[1])                 # ln(y) = ln(a) + ln(x^b)
  output$b [8] <- coef(mod8)[2]                      # ln(y) = ln(a) + b*ln(x)
  output$R2[8] <- rsquare(yg0, output$a[8]*xg0^output$b[8], quiet=quiet)
  } else
  if(!quiet) warning("log/exp/power not fitted (no non-zero positive values in dataset).")
#  9 reciprocal ----------- a/x + b --------------------------------------------
xn0 <- x[x!=0 & y!=0] # xn0: x not zero
yn0 <- y[x!=0 & y!=0]
if(length(xn0)>0 & length(yn0)>0)
  {
  mod9 <- lm( yn0 ~ I(1/xn0) )
  output$a [9] <- coef(mod9)[2]
  output$b [9] <- coef(mod9)[1]
  output$R2[9] <- rsquare(yn0, output$a[9]/xn0 + output$b[9], quiet=quiet)
  } else
  if(!quiet) warning("reciprocal not fitted (no non-zero values in dataset).")
# 10 rational ------------- 1 / (a*x + b) --------------------------------------
y_rational <- 1/y
infin <- which(!is.finite(y_rational))
if(length(infin)!=0)
  {
  if(!quiet) warning("For rational regressions, ",length(infin),
    " zeros were removed from x and y (",round(length(infin)/length(x)*100,1),"%).")
  x_rat <- x[-infin] # x_rat: rational regression values
  y_rat <- y[-infin]
  } else
  {x_rat <- x; y_rat <- y}
if(length(x_rat)>0 & length(y_rat)>0)
  {
  mod10 <- lm( y_rat ~ x_rat)
  output$a [10] <- coef(mod10)[2]
  output$b [10] <- coef(mod10)[1]
  output$R2[10] <- rsquare(y, 1 / (output$a[10]*x + output$b[10]), quiet=quiet)
  }
# 11 exp_4 ---------------- a*e^(b*(x+c))+d ------------------------------------
# 4-parametric exponential distibutions
if(exp_4)
  {
  output_exp4p <- exp4p(x,y, digits=digits)
  if(Poly45) output_exp4p <- cbind(output_exp4p, e=NA,f=NA)
  output[11,] <- output_exp4p[1,] # only include best fit for plotting
  }
#
# 12 hyperbolic ----------- sinh(a*x+b)+c, cosh(), tanh() ----------------------
# yet to add
#
# name output rows -------------------------------------------------------------
rownames(output) <- c("linear", "square", "cubic", "poly4", "poly5",
     "logarithmic", "exponential", "power", "reciprocal", "rational", "exp_4p" )#, "hyperbolic")
output$nr <- 1:11
#
# Formulas of fitted functions -------------------------------------------------
output$Formulas[1] <- paste0(yf," = ", ab1(output$a[1]),"*",xf,      ab(output$b[1]) )
output$Formulas[2] <- paste0(yf," = ", ab1(output$a[2]),"*",xf,"^2", ab(output$b[2]),"*",xf,      ab(output$c[2]) )
output$Formulas[3] <- paste0(yf," = ", ab1(output$a[3]),"*",xf,"^3", ab(output$b[3]),"*",xf,"^2", ab(output$c[3]),"*",xf,      ab(output$d[3]) )
if(Poly45){
output$Formulas[4] <- paste0(yf," = ", ab1(output$a[4]),"*",xf,"^4", ab(output$b[4]),"*",xf,"^3", ab(output$c[4]),"*",xf,"^2", ab(output$d[4]),"*",xf,      ab(output$e[4]) )
output$Formulas[5] <- paste0(yf," = ", ab1(output$a[5]),"*",xf,"^5", ab(output$b[5]),"*",xf,"^4", ab(output$c[5]),"*",xf,"^3", ab(output$d[5]),"*",xf,"^2", ab(output$e[5]),"*",xf, ab(output$f[5]) )
} # end if Poly45
output$Formulas[6] <- paste0(yf," = ", ab1(output$a[6]),"*log10(",xf, ")", ab(output$b[6]) )
output$Formulas[7] <- paste0(yf," = ", ab1(output$a[7]),"*e^(",           ab1(output$b[7]), "*", xf, ")" )
output$Formulas[8] <- paste0(yf," = ", ab1(output$a[8]),"*", xf, "^",     ab1(output$b[8]) )
output$Formulas[9] <- paste0(yf," = ", ab1(output$a[9]),"/",xf,            ab(output$b[9]) )
output$Formulas[10]<- paste0(yf," = 1/( ", ab1(output$a[10]),              ab(output$b[10]), "*", xf, " )" )
if(exp_4)
output$Formulas[11]<- paste0(yf, " = ", ab1(output$a[11]), "*e^(", ab1(output$b[11]), "*(", xf, ab(output$c[11]), "))", ab(output$d[11]))
#
# edit Rsquared columns in output table ----------------------------------------
output$R2full <- output$R2
output$R2 <- round(output$R2, digits)
ord <- order(output$R2full, decreasing=TRUE) # descending order of goodness of fit, for legend
#
# plot data and functions ------------------------------------------------------
if(all(is.na(y)) | all(is.na(x))) plot <- FALSE
if(plot & nbest!=0) {
  if(!add)  plot(x, y, las=las, pch=pch, ylab=ylab, xlab=xlab, xlim=xlim, ylim=ylim, col=pcol, ...)
  # select function types that should be drawn:
  todraw <- rep(FALSE, 12)
  if(missing(R2min) & is.null(selection)) todraw[ord[1:nbest]] <- TRUE
  if(!is.null(selection)) todraw[selection] <- TRUE
  if(!missing(R2min)) todraw[output$R2full>=R2min] <- TRUE
  #
  xdraw <- seqR(par("usr")[1:2], len=200)
  xdrawtab <- data.frame(x=xdraw) # colnames(xdrawtab) <- as.character(xf) # not necessary, as poly uses "x" (the one in the function environment)
  if(is.null(col)) col <- c("black", "red", "green3", "chartreuse", "forestgreen", "blue", "cyan", "magenta", "yellow", "gray", "orange", "deeppink")
  col <- rep(col, length=12)
  #
  if(todraw[1]) lines(xdraw, predict( mod1 , xdrawtab ),          col=col[1], lwd=lwd[1], lty=lty[1]) # 1 linear
  if(todraw[2]) lines(xdraw, predict( mod2 , xdrawtab ),          col=col[2], lwd=lwd[2], lty=lty[2]) # 2 square
  if(todraw[3]) lines(xdraw, predict( mod3 , xdrawtab ),          col=col[3], lwd=lwd[3], lty=lty[3]) # 3 cubic
  if(Poly45){
  if(todraw[4]) lines(xdraw, predict( mod4 , xdrawtab ),          col=col[4], lwd=lwd[4], lty=lty[4]) # 4 polynomial 4th degree
  if(todraw[5]) lines(xdraw, predict( mod5 , xdrawtab ),          col=col[5], lwd=lwd[5], lty=lty[5]) # 5 polynomial 5th degree
  } # end if Poly45
  if(all(xdraw<=0)) if(!quiet) warning("no logarithmic regression could be done, as there are no positive x values")
  xd2 <- xdraw[xdraw>0]
  if(todraw[6]) lines(xd2,   output$a[6]*log10(xd2)+output$b[6],  col=col[6], lwd=lwd[6], lty=lty[6]) # 6 logarithmic
  if(todraw[7]) lines(xdraw, output$a[7]*exp(output$b[7]*xdraw),  col=col[7], lwd=lwd[7], lty=lty[7]) # 7 exponential
  if(todraw[8]) lines(xdraw, output$a[8]*xdraw^output$b[8],       col=col[8], lwd=lwd[8], lty=lty[8]) # 8 power (Potenz)
  if(todraw[9]) lines(xdraw, output$a[9]/xdraw+output$b[9],       col=col[9], lwd=lwd[9], lty=lty[9]) # 9 reciprocal
  if(todraw[10])lines(xdraw, 1/(output$a[10]*xdraw+output$b[10]), col=col[10],lwd=lwd[10],lty=lty[10])# 10 rational
  if(exp_4)
  if(todraw[11])lines(xdraw, output$a[11]*exp(output$b[11]*(xdraw+output$c[11]))+output$d[11], col=col[11],lwd=lwd[11], lty=lty[11]) # 11 exp_4par
  # if(todraw[12])lines(xdraw, output$a[12]*cosh(x)+output$b[12], col=col[12],lwd=lwd[12], lty=lty[12]) # 12 hyperbolic
# prepare and write legend -----------------------------------------------------
if(legend) {
  fForms <- c("a*x + b", "a*x^2 + b*x + c", "a*x^3 + b*x^2 + c*x + d",
     "a*x^4 + b*x^3 + c*x^2 + d*x + e", "a*x^5 + b*x^4 + c*x^3 + d*x^2 + e*x + f",
     "a*log(x) + b", "a*e^(b*x)", "a*x^b", "a/x + b", "1 / (a*x + b)", "a*e^(b*(x+c)) + d")
  #
  legendlabel <- if(legendform=="full") output$Formulas else
                 if(legendform=="form") fForms else
                 if(legendform=="nameform") paste(rownames(output), "  ", fForms) else
                 if(legendform=="name") rownames(output) else
                 stop("wrong legendform. Use 'full', 'form', 'nameform', or 'name'.")
  # remove entries for Poly45 and exp_4 if necessary:
  ###ord <- ord[ !is.na(output$a)]
  if(!Poly45) ord <- ord[-(length(ord)-0:1)] # remove last two (Poly4 and 5) from legend
  if(!exp_4)  ord <- ord[-length(ord)] # remove last one (exp_4) from legend
  ord <- ord[ord %in% which(todraw) ] # keep only the ones that are to be plotted
  #
  legargdefaults <- list(x="top", bty="n", cex=0.8, text.col=col[ord],
       legend=paste(sprintf(paste0("%.",digits,"f"), output$R2), "  ", legendlabel)[ord])
  do.call(graphics::legend, args=owa(legargdefaults, legargs, "legend"))
  } # if legend end
  } # if plot end
#
# final output -----------------------------------------------------------------
if(exp_4) output <- rbind(output, output_exp4p[-1,]) # best fit is already included
if(!exp_4)  output <- output[-11,   ] # remove excess rows
if(!Poly45) output <- output[-(4:5),] # remove excess rows
#
output <- output[order(output$R2full, decreasing=TRUE),]
#
if(ncolumns >10) ncolumns <- 10
if(ncolumns <0) ncolumns <- 0
if(ncolumns >8 & !Poly45) ncolumns <- 8
#
if(ncolumns !=0) return(output[,1:ncolumns])
#
} # Function end ---------------------------------------------------------------
