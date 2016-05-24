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
#       we now prefer "dayrange" to "yearRng" for the weather data.
#    -- three of us wrote the book, and the person preparing these scripts
#       might not be the person who wrote the text
#  Moreover, we expect to augment and modify these command scripts from time
#  to time as we get new data illustrating new things, add functionality
#  to the package, or just for fun.

###
### ch. Chapter 7  Exploring Variation: Functional Principal
###                and Canonical Components analysis
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
## Section 7.1 An Overview of Functional PCA
##
#  (no computations in this section)

##
## Section 7.2 PCA with Function pca.fd
##

# Section 7.2.1 PCA of the Log Precipitation Data

#  Create logprec.fd, copy from chapter 6

logprecav = CanadianWeather$dailyAv[
              dayOfYearShifted, , 'log10precip']
dayrange  = c(0,365)
daybasis  = create.fourier.basis(dayrange, 365)

Lcoef        = c(0,(2*pi/diff(dayrange))^2,0)
harmaccelLfd = vec2Lfd(Lcoef, dayrange)

lambda      = 1e6
fdParobj    = fdPar(daybasis, harmaccelLfd, lambda)
logprec.fit = smooth.basis(day.5, logprecav, fdParobj)
logprec.fd  = logprec.fit$fd

#  do principal component analysis with 2 components

nharm = 2
logprec.pcalist = pca.fd(logprec.fd, nharm)

print(logprec.pcalist$values[1:4])

# Figure 7.1

plot.pca.fd(logprec.pcalist)

#  The expansion supplied by the function is too large,
#  and here we supply a smaller value, 0.5

plot(logprec.pcalist, expand=.5)

# Figure 7.2

logprec.rotpcalist = varmx.pca.fd(logprec.pcalist)

plot.pca.fd(logprec.rotpcalist, expand=.5)

# Figure 7.3

rotpcascores = logprec.rotpcalist$scores

plot(rotpcascores[,1], rotpcascores[,2], type="p", pch="o",
     xlab="Rotated Harmonic I", ylab="Rotated Harmonic II")

# Section 7.2.2 PCA of Log Precipitation Residuals

# logprecres = residuals from
# the smooths of the log precipitation curves in Chapter 5.

logprecmat = eval.fd(day.5, logprec.fd)
logprecres = logprecav - logprecmat

# Figure 7.4

logprecres.fd = smooth.basis(day.5, logprecres,
    fdParobj)$fd
plot(logprecres.fd, lwd=2, col=4, lty=1, cex=1.2,
     xlim=c(0,365), ylim=c(-0.07, 0.07),
     xlab="Day", ylab="Residual (log 10 mm)")

# Figure 7.5

logprec.pca1 = pca.fd(logprecres.fd, 1)
plot(logprec.pca1, expand=0.01)

##
## Section 7.3 More Functional PCA Features
##

#  (no computations in this section)

##
## Section 7.4 PCA of joint X-Y Variation in Handwriting
##

#  Define time values and order 6 spline basis with 105 basis functions
#  This places a knot at every 23rd observation point, and is found to
#  correspond closely to spline smoothing results.

fdabasis = create.bspline.basis(c(0, 2300), 105, 6)
fdatime = seq(0, 2300, len=1401)

#  set up the functional data structure

fdafd = smooth.basis(fdatime, handwrit, fdabasis)$fd
fdafd$fdnames[[1]] = "Milliseconds"
fdafd$fdnames[[2]] = "Replications"
fdafd$fdnames[[3]] = list("X", "Y")

#  plot the data

op <- par(mfrow=c(2,1))
plot(fdafd)
par(op)

#  a principal components analysis

nharm = 3
fdapcaList = pca.fd(fdafd, nharm)

plot.pca.fd(fdapcaList, expand=.2)

fdarotpcaList = varmx.pca.fd(fdapcaList)
plot.pca.fd(fdarotpcaList, expand=.2)

fdaeig = fdapcaList$values
neig = 12
x = matrix(1,neig-nharm,2)
x[,2] = (nharm+1):neig
y = as.matrix(log10(fdaeig[(nharm+1):neig]))
c = lsfit(x,y,int=FALSE)$coef

# Figure 7.6

op <- par(mfrow=c(1,1),cex=1.2)
plot(1:neig, log10(fdaeig[1:neig]), "b",
     xlab="Eigenvalue Number",
     ylab="Log10 Eigenvalue")
lines(1:neig, c[1]+ c[2]*(1:neig), lty=2)
par(op)

# Figure 7.7 varimax rotation

#  set up mean function

fdameanfd  = mean(fdafd)
fdameanmat = eval.fd(fdatime, fdameanfd)

#  evaluate the harmonics

harmfd  = fdarotpcaList$harm
harmmat = eval.fd(fdatime, harmfd)

fdapointtime = seq(0,2300,len=201)
fdameanpoint = eval.fd(fdapointtime, fdameanfd)
harmpointmat = eval.fd(fdapointtime, harmfd)

fac = 0.1
harmplusmat = array(0,c(201,3,2))
harmminsmat = array(0,c(201,3,2))
for (j in 1:3) {
    harmplusmat[,j,] = fdameanpoint[,1,] + fac*harmpointmat[,j,]
    harmminsmat[,j,] = fdameanpoint[,1,] - fac*harmpointmat[,j,]
}

j=3
    plot(fdameanmat[,1,1]-0.035,  fdameanmat[,1,2], "l", lwd=2,
         xlim=c(-0.075,0.075), ylim=c(-0.04, 0.04),
         xlab="", ylab="")
    lines(harmplusmat[,j,1]-0.035, harmplusmat[,j,2], lty=2)
    lines(harmminsmat[,j,1]-0.035, harmminsmat[,j,2], lty=2)
j=2
    lines(fdameanmat[,1,1]+0.035,  fdameanmat[,1,2],  lty=1, lwd=2)
    lines(harmplusmat[,j,1]+0.035, harmplusmat[,j,2], lty=2)
    lines(harmminsmat[,j,1]+0.035, harmminsmat[,j,2], lty=2)


##
## Section 7.5 Exploring Functional Covariation
##             with Canonical Correlation Analysis
##

#  set up temp.fd

tempav = CanadianWeather$dailyAv[
              dayOfYearShifted, , 'Temperature.C']

lambda   = 1e2
fdParobj = fdPar(daybasis, harmaccelLfd, lambda)
temp.fd  = smooth.basis(day.5, tempav, fdParobj)$fd
temp.fd$fdnames = list("Day (July 2 to June 30)",
                       "Weather Station",
                       "Mean temperature (deg. C)")

ccafdPar = fdPar(daybasis, 2, 5e6)
ccalist  = cca.fd(temp.fd, logprec.fd, 3, ccafdPar, ccafdPar)

ccawt.temp    = ccalist$ccawtfd1
ccawt.logprec = ccalist$ccawtfd2
corrs         = ccalist$ccacorr

print(corrs[1:3])
#  [1] 0.9139817 0.6194850 0.3495515

ccawtmat.temp    = eval.fd(day.5, ccawt.temp)
ccawtmat.logprec = eval.fd(day.5, ccawt.logprec)

#  Figure 7.8

plot(day.5, ccawtmat.temp[,1], type='l', lwd=2, cex=2,
     xlab="Day (July 1 to June 30)",
     ylab="Canonical Weight Functions")
lines(day.5, ccawtmat.logprec[,1], lty=2, lwd=2)
lines(dayrange, c(0, 0), lty=3)
legend("bottomleft", c("Temp.", "Log Prec."), lty=c(1,2))

#  Figure 7.9

ccascr.temp    = ccalist$ccavar1
ccascr.logprec = ccalist$ccavar2

placeindex = c(35,30,31,19,33,25,24,17,16,8,14,12,15,10,27,6,1,29)

plot(ccascr.temp[,1], ccascr.logprec[,1], type="p", pch="*", cex=2,
     xlim=c(-40,80),
     xlab="Temperature Canonical Weight",
     ylab="Log Precipitation Canonical Weight")
text(ccascr.temp[placeindex,1]+10, ccascr.logprec[placeindex,1],
     CanadianWeather$place[placeindex])

##
## Section 7.6 Details for the pca.fd and cca.fd Functions
##
help(pca.fd)
help(cca.fd)

##
## Section 7.7 Some Things to Try
##
# (exercises for the reader)

##
## Section 7.8 More to Read
##
