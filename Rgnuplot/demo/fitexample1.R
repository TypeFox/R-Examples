#minpack.lm
###### example 1 from minpack.lm use a seed for reproducibility
set.seed(0)
## values over which to simulate data
x <- seq(0, 5, length = 100)
## model based on a list of parameters
getPred <- function(parS, xx) parS$a * exp(xx * parS$b) + parS$c
## parameter values used to simulate data
pp <- list(a = 9, b = -1, c = 6)
## simulated data, with noise
simDNoisy <- getPred(pp, x) + rnorm(length(x), sd = 0.1)
## plot data
plot(x, simDNoisy, main = "data")
## residual function
residFun <- function(p, observed, xx) observed - getPred(p, xx)
## starting values for parameters
parStart <- list(a = 3, b = -0.001, c = 1)
## perform fit
nls.out <- minpack.lm::nls.lm(par = parStart, fn = residFun, observed = simDNoisy, xx = x, control = minpack.lm::nls.lm.control(nprint = 1))
## plot model evaluated at final parameter estimates
lines(x, getPred(as.list(coef(nls.out)), x), col = 2, lwd = 2)
## summary information on parameter estimates
summary(nls.out)
# Sys.sleep(5);dev.off()

# fit using Rgnuplot
require(Rgnuplot)
# create temporary files for the data and the fit results
tmpfile <- tempfile()
logfile <- tempfile()
write(t(cbind(x, simDNoisy)), tmpfile, ncolumns = 2)
# Initialize the gnuplot handle
h1 <- Gpinit()
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1, "a=9\nb=-1\nc=6\ny(x)=a * exp(b * x) + c\nset fit logfile \"" %s% logfile %s% "\"\nfit y(x) \"" %s% tmpfile %s% "\" u 1:2 via a,b,c\nset xlabel \"x\"\nset ylabel \"y\"\nset nokey\n#set terminal postscript eps color;set output \"fitexample1.eps\"\nplot y(x), \"" %s% 
    tmpfile %s% "\"")
# show the fit progress from the log file, summary information
cat(Gpfile2string(logfile))
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
