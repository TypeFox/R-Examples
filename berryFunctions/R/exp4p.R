#' 4-parametric exponential function
#' 
#' Fits an exponential function of the form a*e^(b*(x+c))+d
#' 
#' @details This is mainly a building block for mReg
#' 
#' @return Data.frame with the 4 parameters for each \code{\link{optim}} method
#' @note Optim can be slow! It refers to the functions rmse and rsquare, also in this package. 
#'       L-BFGS-B needs finite values. In case it doesn't get any 
#'       with the initial parameters (as in the first example Dataset), 
#'       it trys again with the parameters optimized via Nelder Mead.
#' @author Berry Boessenkool, \email{berry-b@@gmx.de}, 2012-2013, outsourced from mReg in July 2014
#' @seealso \code{\link{mReg}}, \code{\link{lm}}
#' @keywords regression nonlinear
#' @export
#' @examples
#' 
#' # exponential decline of temperature of a mug of hot chocolate
#' tfile <- system.file("extdata/Temp.txt", package="berryFunctions")
#' temp <- read.table(tfile, header=TRUE, dec=",")
#' head(temp)
#' plot(temp)
#' temp <- temp[-20,] # missing value - rmse would complain about it
#' x <- temp$Minuten
#' y <- temp$Temp 
#' rm(tfile, temp)
#' 
#' exp4p(x,y, plot=TRUE)
#' # y=49*e^(-0.031*(x - 0  )) + 25 correct, judged from the model:
#' # Temp=T0 - Te *exp(k*t) + Te     with    T0=73.76,  Tend=26.21, k=-0.031
#' # optmethod="Nelder-Mead"  # y=52*e^(-0.031*(x + 3.4)) + 26 wrong
#' 
#' @param x,y x and y Data
#' @param digits significant digits for rounding R^2. DEFAULT: 2
#' @param plot plot data and fitted functions? DEFAULT: FALSE
#' @param las label axis style, see \code{\link{par}}. DEFAULT: 1
#' @param col 6 colors for lines and legend texts. DEFAULT: 1:6
#' @param legarg Arguments passed to \code{\link{legend}}. DEFAULT: NULL
#' @param \dots further graphical parameters passed to \code{\link{plot}}
#' 
exp4p <- function(
    x, y,
    digits=2,
    plot=FALSE,
    las=1,
    col=1:6,
    legarg=NULL,
    ...)
{
# Prepare Output Table
output <- as.data.frame(matrix(NA, ncol=8, nrow=6 ))
colnames(output) <- c("nr","R2","Formulas","R2full", letters[1:4] )
rownames(output) <- c("ini", "N-M", "BFGS", "CG", "SANN", "L--B")
#
# initial parameters via lm of values relocated to first quadrant
init_c <- -min(x, na.rm=TRUE)
init_d <- min(y, na.rm=TRUE)
y_ini <-  y - init_d  + 0.05*abs(init_d)   #; y_ini[y_ini==0] <- 0.001
x_ini <-  x + init_c
init_model <- lm( log(y_ini) ~ x_ini ) # exponential model
init_a <- exp(coef(init_model)[1])
init_b <- coef(init_model)[2]
param <- c(a=init_a, b=init_b, c=init_c, d=init_d)
names(param) <- letters[1:4]
#
# Exponential function to be fitted via optim
expfun <- function(p, x) p[1]*exp(p[2]*(x+p[3]))+p[4]
# function returning one value to be minimized via optim
minfun <- function(p) rmse(y, expfun(p, x=x)) # Root Mean Square Error
#
# Fitting of parameters with different methods
output[1, 5:8] <- param
output[2, 5:8] <- optim(par=param, fn=minfun, method="Nelder-Mead")$par
output[3, 5:8] <- optim(par=param, fn=minfun, method="BFGS")$par
output[4, 5:8] <- optim(par=param, fn=minfun, method="CG")$par
output[5, 5:8] <- optim(par=param, fn=minfun, method="SANN")$par
opt_L      <- try(optim(par=param, fn=minfun, method="L-BFGS-B")$par, silent=TRUE)
for(i in 2:5) if(inherits(opt_L, "try-error"))
  opt_L<- try(optim(output[i,5:8], fn=minfun, method="L-BFGS-B")$par, silent=TRUE)
if(inherits(opt_L, "try-error")) opt_L <- list(par=c(a=NA,b=NA,c=NA,d=NA))
output[6, 5:8] <- opt_L
#
# R squared values
output$R2full <- sapply(1:6, function(i) rsquare(y, expfun(as.numeric(output[i, 5:8]), x=x)))
output$R2 <- round(output$R2full, digits)
# descending order of goodness of fit, for legend
output <- output[ order(output$R2full, decreasing=TRUE) , ]
#
# plot data and function fit ------------------------------------------------------
if(plot)
  {
  plot(x, y, las=las, ...)
  xdraw <- seqR(par("usr")[1:2], len=200)
  for(i in 1:6)
  lines(xdraw, output$a[i]*exp(output$b[i]*(xdraw+output$c[i]))+output$d[i], col=col[i])
  do.call(legend, owa(list(x="topright", legend=rownames(output), col=col, lty=1), legarg, "col"))
}
#
# Return Output
output
} # Function end ---------------------------------------------------------------

