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
### ch. 5.  Smoothing: Computing Curves from Noisy Data
###

#  load the fda package

library(fda)

#  display the data files associated with the fda package

data(package='fda')

#  start the HTML help system if you are connected to the Internet, in
#  order to open the R-Project documentation index page in order to obtain
#  information about R or the fda package.

help.start()

##
## Section 5.1.  Regression Splines: Smoothing by Regression Analysis
##

#  -------------------  Smoothing the growth data  ------------------------

#  define the range of the ages and set up a fine mesh of ages

ageRng  = c(1,18)
age     = growth$age
agefine = seq(1,18,len=501)

#  set up order 6 spline basis with 12 basis functions for
#  fitting the growth data so as to estimate acceleration

nbasis = 12;
norder =  6;
heightbasis12 = create.bspline.basis(ageRng, nbasis, norder)

#  fit the data by least squares

basismat   = eval.basis(age, heightbasis12)
heightmat  = growth$hgtf
heightcoef = lsfit(basismat, heightmat, intercept=FALSE)$coef

#  fit the data using function smooth_basis, which does the same thing.

heightList = smooth.basis(age, heightmat, heightbasis12)
heightfd   = heightList$fd
height.df  = heightList$df
height.gcv = heightList$gcv

heightbasismat = eval.basis(age, heightbasis12)
y2cMap         = solve(crossprod(heightbasismat),
                       t(heightbasismat))

##
## Section 5.2.  Data Smoothing with Roughness Penalties
##

# section 5.2.2 The Roughness Penalty Matrix R

#  ---------------  Smoothing the Canadian weather data  ------------------

#  define a Fourier basis for daily temperature data

yearRng   = c(0,365)
nbasis    = 65
tempbasis = create.fourier.basis(yearRng,nbasis)

#  define the harmonic acceleration operator

#  When the coefficients of the linear differential operator
#  are constant, the Lfd object can be set up more simply
#  by using function vec2Lfd as follows:

#  The first argument is a vector of coefficients for the
#  operator, and the second argument is the range over which
#  the operator is defined.

harmaccelLfd = vec2Lfd(c(0,(2*pi/365)^2,0), yearRng)

#  compute the penalty matrix R

Rmat = eval.penalty(tempbasis, harmaccelLfd)

# section 5.2.4 Defining Smoothing by Functional Parameter Objects

#  -------  Smoothing the growth data with a roughness penalty  -----------

#  set up a basis for the growth data
#  with knots at ages of height measurement

norder      = 6
nbasis      = length(age) + norder - 2
heightbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#  define a functional parameter object for smoothing

heightLfd    = 4
heightlambda = 0.01
heightfdPar  = fdPar(heightbasis, 4, 0.01)

#  smooth the data

heightfdSmooth = smooth.basis(age, heightmat, heightfdPar)
heightfd       = heightfdSmooth$fd

# section 5.2.5 Choosing Smoothing Parameter lambda

loglam         = seq(-6, 0, 0.25)
Gcvsave        = rep(NA, length(loglam))
names(Gcvsave) = loglam
Dfsave         = Gcvsave
for(i in 1:length(loglam)){
  hgtfdPari  = fdPar(heightbasis, Lfdobj=4, 10^loglam[i])
  hgtSm.i    = smooth.basis(age, heightmat, hgtfdPari)
  Gcvsave[i] = sum(hgtSm.i$gcv)
  Dfsave[i]  = hgtSm.i$df
}

# Figure 5.1.

plot(loglam, Gcvsave, 'o', las=1, xlab=expression(log[10](lambda)),
     ylab=expression(GCV(lambda)), lwd=2 )

##
## 5.3.  Case Study: The Log Precipitation Data
##

#  organize data to have winter in the center of the plot

dayOfYearShifted = c(182:365, 1:181)

logprecav = CanadianWeather$dailyAv[
         dayOfYearShifted, , 'log10precip']

#  set up a saturated basis: as many basis functions as observations

nbasis   = 365
daybasis = create.fourier.basis(yearRng, nbasis)

#  set up the harmonic acceleration operator

Lcoef        = c(0,(2*pi/diff(yearRng))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, yearRng)

#  step through values of log(lambda)

loglam        = seq(4,9,0.25)
nlam          = length(loglam)
dfsave        = rep(NA,nlam)
names(dfsave) = loglam
gcvsave       = dfsave
for (ilam in 1:nlam) {
  cat(paste('log10 lambda =',loglam[ilam],'\n'))
  lambda        = 10^loglam[ilam]
  fdParobj      = fdPar(daybasis, harmaccelLfd, lambda)
  smoothlist    = smooth.basis(day.5, logprecav,
                            fdParobj)
  dfsave[ilam]  = smoothlist$df
  gcvsave[ilam] = sum(smoothlist$gcv)
}

# Figure 5.2.

plot(loglam, gcvsave, type='b', lwd=2, ylab='GCV Criterion',
     xlab=expression(log[10](lambda)) )

#  smooth data with minimizing value of lambda

lambda      = 1e6
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprec.fit$fd
fdnames     = list("Day (July 1 to June 30)",
                   "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")
logprec.fd$fdnames = fdnames

#  plot the functional data object

plot(logprec.fd, lwd=2)

# plotfit.fd:  Pauses between plots
# *** --->>> input required (e.g., click on the plot)
#            to advance to the next plot

plotfit.fd(logprecav, day.5, logprec.fd, lwd=2)

##
## Section 5.4 Positive, Monotone, Density
##             and Other Constrained Functions
##

#   ----------------  Positive smoothing of precipitation  ----------------

lambda      = 1e3
WfdParobj   = fdPar(daybasis, harmaccelLfd, lambda)
VanPrec     = CanadianWeather$dailyAv[
  dayOfYearShifted, 'Vancouver', 'Precipitation.mm']
VanPrecPos  = smooth.pos(day.5, VanPrec, WfdParobj)
Wfd         = VanPrecPos$Wfdobj
Wfd$fdnames = list("Day (July 1 to June 30)",
      "Weather Station" = CanadianWeather$place,
                   "Log 10 Precipitation (mm)")

precfit = exp(eval.fd(day.5, Wfd))

plot(day.5, VanPrec, type="p", cex=1.2,
     xlab="Day (July 1 to June 30)",
     ylab="Millimeters",
     main="Vancouver's Precipitation")
lines(day.5, precfit,lwd=2)

#  ------------  5.4.2.1  Monotone smoothing  of the tibia data  -----------

#  set up the data for analysis

day    = infantGrowth[, 'day']
tib    = infantGrowth[, 'tibiaLength']
n      = length(tib)

#  a basis for monotone smoothing

nbasis = 42
Wbasis   = create.bspline.basis(c(1,n), nbasis)

#  the fdPar object for smoothing

Wfd0     = fd(matrix(0,nbasis,1), Wbasis)
WfdPar   = fdPar(Wfd0, 2, 1e-4)

#  smooth the data

result   = smooth.monotone(day, tib, WfdPar)
Wfd      = result$Wfd
beta     = result$beta

#  compute fit and derivatives of fit

dayfine  = seq(1,n,len=151)
tibhat   = beta[1]+beta[2]*eval.monfd(dayfine ,Wfd)
Dtibhat  =        beta[2]*eval.monfd(dayfine, Wfd, 1)
D2tibhat =        beta[2]*eval.monfd(dayfine, Wfd, 2)

#  plot height

op = par(mfrow=c(3,1), mar=c(5,5,3,2), lwd=2)

plot(day, tib, type = "p", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Length (mm)")
lines(dayfine, tibhat, lwd=2)

#  plot velocity

plot(dayfine, Dtibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Velocity (mm/day)")

#  plot acceleration

plot(dayfine, D2tibhat, type = "l", cex=1.2, las=1,
     xlab="Day", ylab='', main="Tibia Acceleration (mm/day/day)")
lines(c(1,n),c(0,0),lty=2)

par(op)

# ---------  5.4.2.2  Monotone smoothing the Berkeley female data  --------

##
##  Compute the monotone smoothing of the Berkeley female growth data.
##

#  set up ages of measurement and an age mesh

age     = growth$age
nage    = length(age)
ageRng  = range(age)
nfine   = 101
agefine = seq(ageRng[1], ageRng[2], length=nfine)

#  the data

hgtf   = growth$hgtf
ncasef = dim(hgtf)[2]

#  an order 6 bspline basis with knots at ages of measurement

norder = 6
nbasis = nage + norder - 2
wbasis = create.bspline.basis(ageRng, nbasis, norder, age)

#  define the roughness penalty for function W

Lfdobj    = 3          #  penalize curvature of acceleration
lambda    = 10^(-0.5)  #  smoothing parameter
cvecf     = matrix(0, nbasis, ncasef)
Wfd0      = fd(cvecf, wbasis)
growfdPar = fdPar(Wfd0, Lfdobj, lambda)

#  monotone smoothing

growthMon = smooth.monotone(age, hgtf, growfdPar)

# (wait for an iterative fit to each of 54 girls)

Wfd        = growthMon$Wfd
betaf      = growthMon$beta
hgtfhatfd  = growthMon$yhatfd

#  Set up functional data objects for the acceleration curves
#  and their mean.  Suffix UN means "unregistered".

accelfdUN     = deriv.fd(hgtfhatfd, 2)
accelmeanfdUN = mean(accelfdUN)

#  plot unregistered curves

par(ask=FALSE)
plot(accelfdUN, xlim=ageRng, ylim=c(-4,3), lty=1, lwd=2,
     cex=2, xlab="Age", ylab="Acceleration (cm/yr/yr)")

# 5.4.3.  Probability density Function

#  -----------  Density function for Regina Precipitation  ----------------

#  plot the empirical quantile function

NR = length(ReginaPrecip)
plot(1:NR, sort(ReginaPrecip), xlab='Rank of rainfall',
     ylab='Ordered daily rainfall (mm)' )

#  set up spline basis for log precipitation density
#  with knots logarithmicaly spaced

sel2.45 = ((2 <= ReginaPrecip) & (ReginaPrecip <= 45))
RegPrec = sort(ReginaPrecip[sel2.45])
N = length(RegPrec)

#  set up spline basis for log precipitation density
#  with knots logarithmicaly spaced

Wknots  = RegPrec[round(N*seq(1/N,1,len=11),0)]
Wnbasis = length(Wknots) + 2
Wbasis  = create.bspline.basis(range(RegPrec),13,4,Wknots)

#  set up the functional parameter object

Wlambda     = 1e-1
WfdPar      = fdPar(Wbasis, 2, Wlambda)

#  estimate the density

densityList = density.fd(RegPrec, WfdPar)
Wfd         = densityList$Wfdobj
C.          = densityList$C

#  plot the density with the knot locations

Zfine = seq(RegPrec[1],RegPrec[N],len=201)
Wfine = eval.fd(Zfine, Wfd)
Pfine = exp(Wfine)/C.

plot(Zfine, Pfine, type='l', lwd=2, xlab='Precipitation (mm)',
     ylab='Probability Density')
abline(v=Wknots, lty='dashed', lwd=2)

##
## Section 5.5 Assessing the Fit to the Log Precipitation Data
##

logprecmat = eval.fd(day.5, logprec.fd)
logprecres = logprecav - logprecmat

#  variance across stations

logprecvar1 = apply(logprecres^2, 1, sum)/35

#  variance across time
logprecvar2 = apply(logprecres^2, 2, sum)/(365-12)

# Figure 5.7

plot(sqrt(logprecvar2), xlab='Station Number',
     ylab='Standard Deviation across Day')
rt = which(CanadianWeather$place %in%
     c("Winnipeg", 'Regina', 'Churchill', 'Montreal', 'St. Johns'))
lft = which(CanadianWeather$place %in%
             c('Yellowknife', 'Resolute', 'Vancouver', 'Iqaluit',
               'Pr. George', 'Pr. Rupert') )
below = which(CanadianWeather$place %in% 'Edmonton')
top = which(CanadianWeather$place %in% 'Halifax')

text(rt, sqrt(logprecvar2[rt]), labels=CanadianWeather$place[rt],
     pos=4)
text(lft, sqrt(logprecvar2[lft]), labels=CanadianWeather$place[lft],
     pos=2)
text(below, sqrt(logprecvar2[below]), labels=CanadianWeather$place[below],
     pos=1)
text(top, sqrt(logprecvar2[top]), labels=CanadianWeather$place[top],
     pos=3)

# Figure 5.8

logstddev.fit = smooth.basis(day.5, log(logprecvar1)/2, fdParobj)
logstddev.fd = logstddev.fit$fd
logprecvar1fit = exp(eval.fd(day.5, logstddev.fd))

plot(day.5, sqrt(logprecvar1), xlab='Day',
     ylab='Standard seviation across stations')
lines(day.5, logprecvar1fit, lwd=2)

##
## Section 5.6 Details for the fdPar Class and smooth.basis Function
##

help(fdPar)
help(smooth.basis)

##
## Section 5.8 Some Things to Try
##
# (exercises for the reader)

##
## Section 5.7 More to Read
##
