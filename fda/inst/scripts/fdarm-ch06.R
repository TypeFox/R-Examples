###
### Ramsay, Hooker & Graves (2009)
### Functional Data Analysis with R and Matlab (Springer)
###

#  Remarks and disclaimers

#  These R commands are either those in this book, or designed to 
#  otherwise illustrate how R can be used in the analysis of functional
#  data.  
#  We do not claim to reproduce the results in the book exactly by these 
#  commands for various reasons, including:
#    -- the analyses used to produce the book may not have been
#       entirely correct, possibly due to coding and accuracy issues
#       in the functions themselves 
#    -- we may have changed our minds about how these analyses should be 
#       done since, and we want to suggest better ways
#    -- the R language changes with each release of the base system, and
#       certainly the functional data analysis functions change as well
#    -- we might choose to offer new analyses from time to time by 
#       augmenting those in the book
#    -- many illustrations in the book were produced using Matlab, which
#       inevitably can imply slightly different results and graphical
#       displays
#    -- we may have changed our minds about variable names.  For example,
#       we now prefer "yearRng" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

###
### ch. 6.  Descriptions of Functional Data
###

library(fda)

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

#  start the HTML help system if you are connected to the Internet, in
#  order to open the R-Project documentation index page in order to obtain
#  information about R or the fda package.

help.start()

##
## Section 6.1 Some Functional Descriptive Statistics
##

#   ----------  Statistics for the log precipitation data  ---------------

#  using 'logprec.fd' computed in fdarm-ch05.R as follows:

#  organize data to have winter in the center of the plot

logprecav = CanadianWeather$dailyAv[
                dayOfYearShifted, , 'log10precip']

#  set up a saturated basis: as many basis functions as observations

yearRng  = c(0,365)
daybasis = create.fourier.basis(yearRng, 365)

#  define the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

#  smooth data with lambda that minimizes GCV

lambda     = 1e6
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fd = smooth.basis(day.5, logprecav, fdParobj)$fd

#  elementary pointwise mean and standard deviation

meanlogprec   = mean(logprec.fd)
stddevlogprec = std.fd(logprec.fd)

# Section 6.1.1 The Bivariate Covariance Function v(s; t)

logprecvar.bifd = var.fd(logprec.fd)

weektime        = seq(0,365,length=53)
logprecvar_mat  = eval.bifd(weektime, weektime,
                              logprecvar.bifd)

# Figure 6.1

persp(weektime, weektime, logprecvar_mat,
      theta=-45, phi=25, r=3, expand = 0.5,
      ticktype='detailed',
      xlab="Day (July 1 to June 30)",
      ylab="Day (July 1 to June 30)",
      zlab="variance(log10 precip)")

contour(weektime, weektime, logprecvar_mat,
        xlab="Day (July 1 to June 30)",
        ylab="Day (July 1 to June 30)")

# Figure 6.2

day5time = seq(0,365,5)
logprec.varmat = eval.bifd(day5time, day5time,
                    logprecvar.bifd)
contour(day5time, day5time, logprec.varmat,
        xlab="Day (July 1 to June 30)",
        ylab="Day (July 1 to June 30)", lwd=2,
        labcex=1)

##
## Section 6.2 The Residual Variance-Covariance Matrix Se
##
#  (no computations in this section)

##
## Section 6.3 Functional Probes rho[xi]
##

# see section 6.5 below

##
## Section 6.4 Phase-plane Plots of Periodic Effects
##

#  -----------  Phase-plane plots for the nondurable goods data  ---------

#  set up a basis for smoothing log nondurable goods index

goodsbasis  = create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)

#  smooth the data using function smooth.basisPar, which
#  does not require setting up a functional parameter object

LfdobjNonDur= int2Lfd(4)

logNondurSm = smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)

# Fig. 6.3 The log nondurable goods index for 1964 to 1967

sel64.67 = ((1964<=index(nondurables)) &
             (index(nondurables)<=1967) )
plot(index(nondurables)[sel64.67],
     log10(nondurables[sel64.67]), xlab='Year',
     ylab='Log10 Nondurable Goods Index', las=1)
abline(v=1965:1966, lty='dashed')

t64.67 = seq(1964, 1967, len=601)
lines(t64.67, predict(logNondurSm, t64.67))

# Section 6.4.1.  Phase-plane Plots Show Energy Transfer
# Figure 6.4.  Phase-plane plot for a simple harmonic function

sin.   = expression(sin(2*pi*x))
D.sin  = D(sin.,  "x")
D2.sin = D(D.sin, "x")

with(data.frame(x=seq(0, 1, length=46)),
     plot(eval(D.sin), eval(D2.sin), type="l",
          xlim=c(-10, 10), ylim=c(-50, 50),
          xlab="Velocity", ylab="Acceleration"), las=1 )
pi.2  = (2*pi)
#lines(x=c(-pi.2,pi.2), y=c(0,0), lty=3)
pi.2.2= pi.2^2
lines(x=c(0,0), y=c(-pi.2.2, pi.2.2), lty="dashed")
lines(x=c(-pi.2, pi.2), y=c(0,0), lty="dashed")
text(c(0,0), c(-47, 47), rep("Max. potential energy", 2))
text(c(-8.5,8.5), c(0,0), rep("Max\nkinetic\nenergy", 2))

# Section 6.4.2 The Nondurable Goods Cycles

# Figure 6.5

phaseplanePlot(1964, logNondurSm$fd)

# sec. 6.4.3.  Phase-Plane Plotting the Growth of Girls

#  -------------  Phase-plane diagrams for the growth data  ---------------

gr.basis = create.bspline.basis(norder=6, breaks=growth$age)
children = 1:10
ncasef   = length(children)
cvecf           = matrix(0, gr.basis$nbasis, ncasef)
dimnames(cvecf) = list(gr.basis$names,
              dimnames(growth$hgtf)[[2]][children])

gr.fd0      = fd(cvecf, gr.basis)
gr.fdPar1.5 = fdPar(gr.fd0, Lfdobj=3, lambda=10^(-1.5))
hgtfmonfd   = with(growth, smooth.monotone(age, hgtf[,children],
                                           gr.fdPar1.5) )
agefine = seq(1,18,len=101)
(i11.7  = which(abs(agefine-11.7) == min(abs(agefine-11.7)))[1])

velffine = predict(hgtfmonfd, agefine, 1);
accffine = predict(hgtfmonfd, agefine, 2);

plot(velffine, accffine, type='n', xlim=c(0, 12), ylim=c(-5, 2),
     xlab='Velocity (cm/yr)', ylab=expression(Acceleration (cm/yr^2)),
     las=1)
for(i in 1:10){
  lines(velffine[, i], accffine[, i])
  points(velffine[i11.7, i], accffine[i11.7, i])
}
abline(h=0, lty='dotted')

##
## Section 6.5 Confidence Intervals for Curves and their Derivatives
##

# sec. 6.5.1.  Two Linear mappings Defining a Probe Value

#  -----------------  Probe values for Canadian weather data  ------------

#  define a probe to emphasize mid-winter

dayvec  = seq(0,365,len=101)
xivec   = exp(20*cos(2*pi*(dayvec-197)/365))
xibasis = create.bspline.basis(c(0,365),13)
xifd    = smooth.basis(dayvec, xivec, xibasis)$fd

plot(xifd)

#  define bases for temperature and precipitation

tempbasis = create.fourier.basis(c(0,365),65)
precbasis = create.fourier.basis(c(0,365),365)

#  probe values for basis functions with respect to xifd

tempLmat = inprod(tempbasis, xifd)
precLmat = inprod(precbasis, xifd)

# sec. 6.5.3.  Confidence Limits for Prince Rupert's Log Precipitation

logprecav = CanadianWeather$dailyAv[
         dayOfYearShifted, , 'log10precip']

# as in section 5.3

#  smooth data with lambda that minimizes GCV getting
#  all of the output up to matrix y2cMap

yearRng      = c(0,365)
daybasis     = create.fourier.basis(yearRng, 365)
Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

lambda     = 1e6
fdParobj   = fdPar(daybasis, harmaccelLfd, lambda)

logprecList = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprecList$fd
fdnames    = list("Day (July 1 to June 30)",
               "Weather Station" = CanadianWeather$place,
               "Log10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  compute the residual matrix and variance vector

logprecmat  = eval.fd(day.5, logprec.fd)
logprecres  = logprecav - logprecmat
logprecvar  = apply(logprecres^2, 1, sum)/(35-1)

#  smooth log variance vector

lambda      = 1e8
resfdParobj = fdPar(daybasis, harmaccelLfd, lambda)
logvar.fd   = smooth.basis(day.5, log(logprecvar), resfdParobj)$fd

#  evaluate the exponentiated log variance vector and
#  set up diagonal error variance matrix SigmaE

varvec      = exp(eval.fd(day.5, logvar.fd))
SigmaE      = diag(as.vector(varvec))

#  compute variance covariance matrix for fit

y2cMap        = logprecList$y2cMap
c2rMap        = eval.basis(day.5, daybasis)
Sigmayhat     = c2rMap %*% y2cMap %*% SigmaE %*%
                t(y2cMap) %*% t(c2rMap)

#  extract standard error function for yhat

logprec.stderr= sqrt(diag(Sigmayhat))

#  plot Figure 6.6

logprec29 = eval.fd(day.5, logprec.fd[29])

plot(logprec.fd[29], lwd=2, ylim=c(0.2, 1.3))
lines(day.5, logprec29 + 2*logprec.stderr,
        lty=2, lwd=2)
lines(day.5, logprec29 - 2*logprec.stderr,
        lty=2, lwd=2)
points(day.5, logprecav[,29])

##
## Section 6.6 Some Things to Try
##
# (exercises for the reader)
