#minpack.lm
###### example 2 from minpack.lm

set.seed(0)  #use a seed for reproducibility
## function to simulate data
f <- function(TT, tau, N0, a, f0) {
    expr <- expression(N0 * exp(-TT/tau) * (1 + a * cos(f0 * TT)))
    eval(expr)
}
## helper function for an analytical gradient
j <- function(TT, tau, N0, a, f0) {
    expr <- expression(N0 * exp(-TT/tau) * (1 + a * cos(f0 * TT)))
    c(eval(D(expr, "tau")), eval(D(expr, "N0")), eval(D(expr, "a")), eval(D(expr, "f0")))
}
## values over which to simulate data
TT <- seq(0, 8, length = 501)
## parameter values underlying simulated data
p <- c(tau = 2.2, N0 = 1000, a = 0.25, f0 = 8)
## get data
Ndet <- do.call("f", c(list(TT = TT), as.list(p)))
## with noise
N <- Ndet + rnorm(length(Ndet), mean = Ndet, sd = 0.01 * max(Ndet))
## plot the data to fit
par(mfrow = c(2, 1), mar = c(3, 5, 2, 1))
plot(TT, N, bg = "black", cex = 0.5, main = "data")
## define a residual function
fcn <- function(p, TT, N, fcall, jcall) (N - do.call("fcall", c(list(TT = TT), as.list(p))))
## define analytical expression for the gradient
fcn.jac <- function(p, TT, N, fcall, jcall) -do.call("jcall", c(list(TT = TT), as.list(p)))
## starting values
guess <- c(tau = 2.2, N0 = 1500, a = 0.25, f0 = 10)
## to use an analytical expression for the gradient found in fcn.jac uncomment jac = fcn.jac
out <- minpack.lm::nls.lm(par = guess, fn = fcn, jac = fcn.jac, fcall = f, jcall = j, TT = TT, N = N, control = minpack.lm::nls.lm.control(nprint = 1))
## get the fitted values
N1 <- do.call("f", c(list(TT = TT), out$par))
## add a blue line representing the fitting values to the plot of data
lines(TT, N1, col = "blue", lwd = 2)
## add a plot of the log residual sum of squares as it is made to decrease each iteration; note that the RSS at the starting parameter values is also stored
plot(1:(out$niter + 1), log(out$rsstrace), type = "b", main = "log residual sum of squares vs. iteration number", xlab = "iteration", ylab = "log residual sum of squares", pch = 21, 
    bg = 2)
## get information regarding standard errors
summary(out)
# end of example 2 from minpack.lm

# fit using Rgnuplot
require(Rgnuplot)
# create temporary files for the data, the fit progress and the fit results
tmpfile <- tempfile()
logfile <- tempfile()
logstderr <- "minpackfit-2.log"
if (file.exists(logstderr)) file.remove(logstderr)
write(t(cbind(TT, N)), tmpfile, ncolumns = 2)
# Initialize the gnuplot handle
h1 <- GpinitSaveStderr(logstderr)
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1, "tau = 2.2\nN0 = 1500\na = 0.25\nf0 = 10\nf(TT) = N0*exp(-TT/tau)*(1 + a*cos(f0*TT))\nset fit logfile \"" %s% logfile %s% "\"\nset dummy TT\nfit f(TT) \"" %s% tmpfile %s% 
    "\" via tau, N0, a, f0")
# wait until the size of the log file is stable
logsize <- -1
while (file.info(logstderr)$size != logsize) {
    Sys.sleep(2)
    logsize <- file.info(logstderr)$size
}
GpfitAllprogress(logstderr)  #shows the progress of the fit algorithm
WSSR <- GpfitProgress(logstderr, "WSSR")  #get all WSSR values
WSSR <- WSSR[-1]  #remove the WSSR from iteration 0
GpMatrixr2gnu(cbind(1:length(WSSR), log(WSSR)), "logWSSR.dat")  #save log(WSSR)

Gpcmd(h1, "#set terminal postscript eps color;set output \"fitexample2.eps\"\nset multiplot layout 2,1\nset title \"Data\"\nset xlabel \"TT\"\nset ylabel \"N\"\nset nokey\nplot \"" %s% 
    tmpfile %s% "\" ls 6 lc rgb \"green\", f(TT) lc rgb \"blue\"\nset title \"log residual sum of squares vs. iteration number\"\nset xlabel \"iteration number\"\nset ylabel \"log residual sum of squares\"\nset nokey\nplot \"logWSSR.dat\" w lp pt 7\nunset multiplot")

print(GpGetvariable(h1, "FIT_NDF"))
print(GpGetvariable(h1, "FIT_WSSR"))
print(GpGetvariable(h1, "FIT_STDFIT"))
# show the fit progress from the log file, summary information
cat(Gpfile2string(logfile))
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
