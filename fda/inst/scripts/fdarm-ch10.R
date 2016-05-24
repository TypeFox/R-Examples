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
### ch.  10.  Linear Models for Functional Responses
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
## Section 10.1  Functional Responses and an Analysis of Variance Model
##

#  ----------------- Climate zone effects for temperature -----------------

#
#  Section 10.1.1 Climate Region Effects on Temperature
#

#  set up the data for the analysis

regions.         = unique(CanadianWeather$region)
p                = length(regions.) + 1
regionList       = vector("list", p)
names(regionList)= c('Canada', regions.)
regionList[[1]]  = c(rep(1,35),0)
for (j in 2:p) {
  xj             = (CanadianWeather$region == regions.[j-1])
  regionList[[j]]= c(xj,1)
}

# tempfd from chapter 9

Lcoef        = c(0,(2*pi/365)^2,0)
harmaccelLfd = vec2Lfd(Lcoef, c(0,365))

tempbasis    = create.fourier.basis(c(0, 365), 65)

lambda       = 1e6
tempfdPar65  = fdPar(tempbasis, harmaccelLfd, lambda)

dayOfYearShifted = c(182:365, 1:181)

tempShifted  = daily$tempav[dayOfYearShifted, ]

tempSmooth65 = smooth.basis(day.5, tempShifted, tempfdPar65)
tempfd       = tempSmooth65$fd

#  augment tempfd by adding a 36th observation with temp(t) = 0

coef    = tempfd$coef
coef36  = cbind(coef,matrix(0,65,1))
temp36fd= fd(coef36,tempbasis,tempfd$fdnames)

#  set up the regression coefficient list

betabasis      = create.fourier.basis(c(0, 365), 11)
betafdPar      = fdPar(betabasis)
betaList       = vector("list",p)
names(betaList)= regions.
for (j in 1:p) betaList[[j]] = betafdPar

#  carry out the functional analysis of variance

fRegressList= fRegress(temp36fd, regionList, betaList)

#  extract the estimated regression coefficients and y-values

betaestList = fRegressList$betaestlist
regionFit   = fRegressList$yhatfd
regions     = c("Canada", regions.)

# Figure 10.1

op          = par(mfrow=c(2,3),cex=1)
for (j in 1:p) plot(betaestList[[j]]$fd, lwd=2,
                    xlab="Day (July 1 to June 30)",
                    ylab="", main=regions[j])
plot(regionFit, lwd=2, col=1, lty=1,
     xlab="Day (July 1 to June 30)", ylab="", main="Prediction")
par(op)

# ------  10.1.2 Trends in Sea Bird Populations on Kodiak Island  ---------

#  select only the data for sites Uyak and Uganik, which have data
#  from 1986 to 2005, except for 1998

sites = c('Uganik', 'Uyak')
sel   = seabird$Bay %in% sites
UU    = seabird[sel,]

# Drop 2 species with many NAs

NAs       = sapply(UU, function(x)sum(is.na(x)))
NAs.      = which(NAs > 2)
birdindex = (1:15)[-NAs.]
birds     = names(UU)[birdindex]

#  Compute mean counts taken over both sites and transects

meanCounts = matrix(NA, 20, 13)
dimnames(meanCounts) = list(1986:2005, birds)

for(i in 1:20){
  sel = (UU$Year == rownames(meanCounts)[i])
  meanCounts[i, ] = sapply(UU[sel, birds], mean, na.rm=TRUE)
}

selYear   = !is.na(meanCounts[,1])
logCounts = log10(meanCounts[selYear,])

#  time vectors in years and in indices in 1:20

yearObs  = as.numeric(rownames(logCounts))
yearCode = (1:20)[selYear]

shellfishindex = c(1,2,5,6,12,13)
fishindex      = (1:13)[-shellfishindex]

# Figure 10.2

ylim = range(logCounts)
op = par(mfrow=c(2,1), mar=c(2, 4, 4, 1)+.1)

matplot(yearObs, logCounts[, shellfishindex], xlab='', ylab='',
        ylim=ylim, main='Shellfish Diet', type='b', col=1)
meanShellfish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanShellfish[!is.na(meanShellfish)]), lwd=3)
abline(h=0, lty='dotted')

matplot(yearObs, logCounts[, fishindex], xlab='', ylab='',
        ylim=ylim, main='Fish Diet', type='b', col=1)
meanFish = apply(meanCounts[, shellfishindex], 1, mean)
lines(yearObs, log10(meanFish[!is.na(meanFish)]), lwd=3)
abline(h=0, lty='dotted')

par(op)

#  Compute mean counts taken over transects only within sites
#  so we have 2 observations for each bird species each year.
#  Two of these counts are zero, and are replaced by 1/(2*n)

meanCounts2 = matrix(NA, 20, 26)

for(i in 1:20) for (j in 1:2) {
  sel = (UU$Year == rownames(meanCounts)[i] & as.character(UU$Bay) == sites[j])
  meanCountsij = sapply(UU[sel, birds], mean, na.rm=TRUE)
  n = sum(sel)
  if (n > 0) {
    meanCountsij[meanCountsij == 0] = 1/(2*n)
  }
  meanCounts2[i,(j-1)*13+(1:13)] = meanCountsij
}

selYear2   = !is.na(meanCounts2[, 1])
yearCode  = (1:20)[selYear2]
all.equal(yearCode, c(1:12, 14:20))

logCounts2 = log10(meanCounts2[selYear2,])

#  Represent log mean counts exactly with a polygonal basis

birdbasis = create.polygonal.basis(yearCode)
birdlist2 = smooth.basis(yearCode, logCounts2, birdbasis)

birdfd2 = birdlist2$fd

#  -----------------------------------------------------------------
#  After some preliminary analyses we determined that there was no
#  contribution from either site or food*site interaction.
#  Now we use a reduced model with only a feed effect,
#  but we add bird effects, which were seen in the plot to be
#  strong.  Birds are nested within feed groups, and either their
#  effects must sum to zero within each group, or we must designate
#  a bird in each group as a baseline, and provide dummy variables
#  for the remainder.  We opt for the latter strategy.
#  -----------------------------------------------------------------

#  The design matrix contains an intercept dummy variable, a
#  feed dummy variable, and dummy variables for birds, excluding
#  the second bird in each group, which turns out to be the each
#  group's most abundant species, and which is designated as the
#  baseline bird for that group.

# 15 columns for the intercept + diet + 13 bird species
# 26 rows for the 26 (species - bay) combinations

Zmat0 = matrix(0,26,15)

#  Intercept or baseline effect

Intercept = rep(1,26)

#  Crustacean/Mollusc feeding effect:  a contrast between the two groups

foodindex = c(1,2,5,6,12,13)
fooddummy = c(2*rep(1:13 %in% foodindex, 2)-1)

#  Bird effect, one for each species

birddummy = diag(rep(1,13))
birdvarbl = rbind(birddummy,birddummy)

#  fill the columns of the design matrix

Zmat0[,1]    = Intercept
Zmat0[,2]    = fooddummy
Zmat0[,3:15] = birdvarbl

#  Two extra dummy observations are added to the functional data
#  object for log counts, and two additional rows are added to
#  the design matrix to force the bird effects within each diet
#  group to equal 0.

birdfd3 = birdfd2
birdfd3$coefs = cbind(birdfd3$coefs, matrix(0,19,2))

Zmat = rbind(Zmat0, matrix(0,2,15))
Zmat[27,shellfishindex+2] = 1
Zmat[28,     fishindex+2] = 1

p = 15
xfdlist = vector("list",p)
names(xfdlist) = c("const", "diet", birds)
betalist = xfdlist
for (j in 1:p) xfdlist[[j]] = Zmat[,j]

#  set up the functional parameter object for (the regression fns.
#  use cubic b-spline basis for intercept and food coefficients

betabasis1 = create.bspline.basis(c(1,20),21,4,yearCode)
Lfdobj1    = int2Lfd(2);
Rmat1      = eval.penalty(betabasis1, Lfdobj1)
lambda1    = 10
betafdPar1 = fdPar(betabasis1,Lfdobj1,lambda1,TRUE,Rmat1)
betalist[[1]] = betafdPar1
betalist[[2]] = betafdPar1
betabasis2 = create.constant.basis(c(1,20))
betafdPar2 = fdPar(betabasis2)
for (j in 3:15) betalist[[j]] = betafdPar2

birdRegress = fRegress(birdfd3, xfdlist, betalist)
betaestlist = birdRegress$betaestlist

# Figure 10.3 is produced in Section 10.2.2 below
# after estimating the smoothing parameter in Section 10.1.3
#
# Here we plot the regression parameters
# without the confidence intervals.

op = par(mfrow=c(2,1))
plot(betaestlist$const$fd)
plot(betaestlist$diet$fd)
par(op)

##
## Section 10.1.3 Choosing Smoothing Parameters
##

#  Choose the level of smoothing by minimizing cross-validated
#  error sums of squares.

loglam = seq(-2,4,0.25)
SSE.CV = rep(0,length(loglam))
betafdPari = betafdPar1
for(i in 1:length(loglam)){
    print(loglam[i])
    betafdPari$lambda = 10^loglam[i]
    betalisti = betalist
    for (j in 1:2) betalisti[[j]] = betafdPari
    CVi = fRegress.CV(birdfd3, xfdlist, betalisti, CVobs=1:26)
    SSE.CV[i] = CVi$SSE.CV
}

#  Figure 10.4

plot(loglam,SSE.CV,type='b',cex.lab=1.5,cex.axis=1.5,lwd=2,
  xlab='log smoothing parameter',ylab='cross validated sum of squares')

#  Cross-validation is minimized at something like lambda = sqrt(10),
#  although the discontinous nature of the CV function is disquieting.

betafdPar1$lambda = 10^0.5
for (j in 1:2) betalist[[j]] = betafdPar1

y = birdfd3
wt = NULL
CVobs = 1:26
returnMatrix=FALSE

#  carry out the functional regression analysis

fitShellfish.5 = fRegress(birdfd3, xfdlist, betalist)

#  plot regression functions

betanames = list("Intercept", "Food Effect")

birdBetaestlist = fitShellfish.5$betaestlist

op = par(mfrow=c(2,1), cex=1.2)
for (j in 1:2) {
    betaestParfdj = birdBetaestlist[[j]]
    betaestfdj    = betaestParfdj$fd
    betaestvecj   = eval.fd(yearCode, betaestfdj)
	  plot(yearObs, betaestvecj, type="l", lwd=4, col=4,
           xlab="Year", ylab="Temp.",
           main=betanames[[j]])
}
par(op)

#  plot predicted functions

birdYhatfdobj = fitShellfish.5$yhatfdobj

plotfit.fd(logCounts2, yearCode, birdYhatfdobj$fd[1:26])
# *** Click on the plot to advance to the next ...

##
## Section 10.2 Functional Responses with Functional Predictors:
##              The Concurrent Model
##

#  Section 10.2.2 Confidence Intervals for Regression Functions

birdYhatmat = eval.fd(yearCode, birdYhatfdobj$fd[1:26])
rmatb   = logCounts2 - birdYhatmat
SigmaEb = var(t(rmatb))

y2cMap.bird = birdlist2$y2cMap

birdStderrList = fRegress.stderr(fitShellfish.5, y2cMap.bird,
                             SigmaEb)
birdBeta.sdList = birdStderrList$betastderrlist

op = par(mfrow=c(2,1))
plotbeta(birdBetaestlist[1:2], birdBeta.sdList[1:2])
par(op)

# Section 10.2.3 Knee Angle Predicted from Hip Angle

gaittime = seq(0.5,19.5,1)
gaitrange = c(0,20)
gaitfine = seq(0,20,len=101)

harmaccelLfd20 = vec2Lfd(c(0, (2*pi/20)^2, 0), rangeval=gaitrange)
gaitbasis = create.fourier.basis(gaitrange, nbasis=21)

gaitLoglam = seq(-4,0,0.25)
nglam   = length(gaitLoglam)

# First select smoothing for the raw data

gaitSmoothStats = array(NA, dim=c(nglam, 3),
      dimnames=list(gaitLoglam, c("log10.lambda", "df", "gcv") ) )
gaitSmoothStats[, 1] = gaitLoglam

#  loop through smoothing parameters

for (ilam in 1:nglam) {
  gaitSmooth = smooth.basisPar(gaittime, gait, gaitbasis,
                   Lfdobj=harmaccelLfd20, lambda=10^gaitLoglam[ilam])
  gaitSmoothStats[ilam, "df"]  = gaitSmooth$df
  gaitSmoothStats[ilam, "gcv"] = sum(gaitSmooth$gcv)
  # note: gcv is a matrix in this case
}

#  display and plot GCV criterion and degrees of freedom

gaitSmoothStats
plot(gaitSmoothStats[, c(1, 3)], type='b')

#  set up plotting arrangements for one and two panel displays
#  allowing for larger fonts

op = par(mfrow=c(2,1))
par(op)
plot(gaitSmoothStats[, c(1, 3)], type="b", log="y")
plot(gaitSmoothStats[, 1:2], type="b", log="y")

#    GCV is minimized with lambda = 10^(-1.5).

gaitSmooth = smooth.basisPar(gaittime, gait,
       gaitbasis, Lfdobj=harmaccelLfd20, lambda=10^(-1.5))
gaitfd = gaitSmooth$fd

names(gaitfd$fdnames) = c("Normalized time", "Child", "Angle")
gaitfd$fdnames[[3]] = c("Hip", "Knee")

hipfd  = gaitfd[,1]
kneefd = gaitfd[,2]

# Figure 10.5

kneefdMean = mean(kneefd)

op = par(mfrow=c(3,1))
plot(kneefdMean, xlab='', ylab='', ylim=c(0, 80),
     main='Mean Knee Angle', lwd=2)
abline(v=c(7.5, 14.7), lty='dashed')
plot(deriv(kneefdMean), xlab='', ylab='',
     main='Knee Angle Velocity', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')
plot(deriv(kneefdMean, 2), xlab='', ylab='',
     main='Knee Angle Acceleration', lwd=2)
abline(v=c(7.5, 14.7), h=0, lty='dashed')
par(op)

# Figure 10.6

phaseplanePlot(gaitfine, kneefdMean,
               labels=list(evalarg=gaittime, labels=1:20),
               xlab='Knee Velocity', ylab='Knee Acceleration')

# Set up a  functional linear regression

xfdlist   = list(const=rep(1,39), hip=hipfd)
betafdPar = fdPar(gaitbasis, harmaccelLfd20)
betalist  = list(const=betafdPar, hip=betafdPar)

gaitRegress= fRegress(kneefd, xfdlist, betalist)

# Figure 10.7
op = par(mfrow=c(2,1))

# Intercept

betaestlist = gaitRegress$betaestlist
kneeIntercept = predict(betaestlist$const$fd, gaitfine)

# mean knee angle

kneeMean = predict(kneefdMean, gaitfine)

# Plot intercept & mean knee angle
ylim1 = range(kneeIntercept, kneeMean)
plot(gaitfine, kneeIntercept, ylim=ylim1, lwd=2,
     main="Intercept and Mean Knee Angle", type='l',
     xlab='', ylab='')
lines(gaitfine, kneeMean, lty='dashed')
abline(h=0, v=c(7.5, 14.7), lty='dashed')

# Hip coefficient

hipCoef = predict(betaestlist$hip$fd, gaitfine)

# Squared multiple correlation

kneehatfd = gaitRegress$yhatfd$fd
kneehatmat = eval.fd(gaittime, kneehatfd)
resmat. = gait[,,'Knee Angle'] - kneehatmat
SigmaE = cov(t(resmat.))

kneefinemat   = eval.fd(gaitfine, kneefd)
kneemeanvec   = eval.fd(gaitfine, mean(kneefd))
kneehatfinemat= eval.fd(gaitfine, kneehatfd)
resmat        = kneefinemat - kneehatfinemat
ncurve        = dim(gait)[2]
resmat0 = kneefinemat - kneemeanvec %*% matrix(1,1,ncurve)
SSE0 = apply((resmat0)^2, 1, sum)
SSE1 = apply(resmat^2, 1, sum)
knee.R2 = (SSE0-SSE1)/SSE0

# Plot Hip Coefficient & Squared Multiple Correlation

ylim2=c(0, max(hipCoef, knee.R2))
plot(gaitfine, hipCoef, lwd=2, xlab='', ylab='', ylim=ylim2, type='l',
     main='Hip Coefficient and Squared Multiple Correlation')
abline(v=c(7.5, 14.7), lty='dashed')
lines(gaitfine, knee.R2, lty='dashed')

# Figure 10.8

gaitbasismat = eval.basis(gaitfine, gaitbasis)
y2cMap = gaitSmooth$y2cMap

fRegressList1 = fRegress(kneefd, xfdlist, betalist,
                         y2cMap=y2cMap, SigmaE=SigmaE)

fRegressList2 = fRegress.stderr(fRegressList1, y2cMap, SigmaE)
betastderrlist = fRegressList2$betastderrlist

op = par(mfrow=c(2,1))
plotbeta(betaestlist, betastderrlist, gaitfine)
par(op)

# Figure 10.9

xfdlist2 = list(const=rep(1,39), hip=deriv(hipfd, 2))
kneefd.accel = deriv(kneefd, 2)
gaitAccelRegr = fRegress(kneefd.accel, xfdlist2, betalist)

gaitt3 = seq(0, 20, length=401)
beta.hipFine = predict(gaitAccelRegr$betaestlist$hip$fd, gaitt3)

plot(gaitt3, beta.hipFine, type ='l', ylim=c(0, max(beta.hipFine)),
     xlab='', ylab='Hip acceleration and squared multiple correlation',
     lwd=2)
abline(v=c(7.5, 14.7), lty='dashed')

# Squared multiple correlation

kneeAccel.pred = predict(gaitAccelRegr$yhatfd$fd, gaitt3)
kneeAccel.     = predict(kneefd.accel, gaitt3)

#MS.pred0 = sd(t(kneeAccel.pred))^2
# warning: sd(matrix) depricated;  use apply(*.2, sd)
MS.pred = apply(kneeAccel.pred, 1, var)
#MS.accelfd0 = sd(t(kneeAccel.))^2
MS.accelfd = apply(kneeAccel., 1, var)
kneeAccel.R2 = (MS.pred / MS.accelfd)

lines(gaitt3, kneeAccel.R2, lty='dashed', lwd=2)

##
## Section 10.3 Beyond the Concurrent Model
##
#  (no computations in this section)

##
## Section 10.4 A Functional Linear Model for Swedish Mortality
##

#  The Swedish mortality data are proprietary, so these data are not
#  included in the fda package.  To use the following code, you must
#  first obtain the data.  The following details are provided to help
#  you to obtain these data.

#  Mortality data for 37 countries can be obtained from the
#  Human Mortality Database (http://www.mortality.org/).
#  For example, the Swedish mortality data can be found at
#  (http://www.mortality.org/cgi-bin/hmd/country.php?cntr=SWE&level=1).

#  Citation:

#  Human Mortality Database. University of California, Berkeley (USA),
#  and Max Planck Institute for Demographic Research (Germany).

#   Two data objects are required for these analyses:

#  SwedeMat:  a dataframe object with 81 rows and 144 columns
#             containing the log hazard values for ages 0 through 80
#             and years 1751 through 1884
#  Swede1920: a vector object containing log hazard values for 1914


## The function readHMD will download and reformat entries in the databases
# lifetables for you which we can then use. Note that this function requires
# the packages RCurl which has sometimes been unavailable on Windows computers.
# If this is the case, the desired entries can be found in the third column
# ('qx') of the Swedish Female Lifetable (file fltcoh_1x1.txt) which must then
# be reformatted to a matrix giving age (rows) and year of birth (columns).

if(FALSE){

# Enter here your USERNAME and PASSWORD
# for www.mortality.org
  USERNAME <- 'jillUser'
  PASSWORD <- 'JUpw123!'

  Sweden = readHMD(USERNAME,PASSWORD,'SWE',ltCol='q')

  SwedeMat = log(Sweden$y[1:81,])

  Swede1920 = SwedeMat[,164]
  SwedeLogHazard = SwedeMat[,1:44]

# SwedeLogHazard = as.matrix(SwedeMat)

  dimnames(SwedeLogHazard)[[2]] <- paste('b', 1751:1894, sep='')

# Figure 10.10

  Fig10.10data = cbind(SwedeLogHazard[, c('b1751', 'b1810', 'b1860')],
                       Swede1920)

  SwedeTime = 0:80;
  SwedeRng = c(0,80);

  matplot(SwedeTime, Fig10.10data,
          type='l',lwd=2,xlab='age',ylab='log Hazard',col=1,
          cex.lab=1.5,cex.axis=1.5)

#  smooth the log hazard observations

  nbasis = 85
  norder = 6
  SwedeBasis = create.bspline.basis(SwedeRng, nbasis, norder)

  D2fdPar = fdPar(SwedeBasis, lambda=1e-7)

  SwedeLogHazfd = smooth.basis(SwedeTime, SwedeLogHazard, D2fdPar)$fd

# The following requires manually clicking on the plot
# for each of 144 birth year cohorts

  plotfit.fd(SwedeLogHazard,SwedeTime,SwedeLogHazfd)

# Set up for the list of regression coefficient fdPar objects

  nbasis     = 23
  SwedeRng   = c(0,80)
  SwedeBetaBasis = create.bspline.basis(SwedeRng,nbasis)

  SwedeBeta0Par = fdPar(SwedeBetaBasis, 2, 1e-5)

  SwedeBeta1fd  = bifd(matrix(0,23,23), SwedeBetaBasis, SwedeBetaBasis)

  SwedeBeta1Par = bifdPar(SwedeBeta1fd, 2, 2, 1e3, 1e3)

  SwedeBetaList = list(SwedeBeta0Par, SwedeBeta1Par)

#  Define the dependent and independent variable objects

  NextYear = SwedeLogHazfd[2:144]
  LastYear = SwedeLogHazfd[1:143]

#  Do the regression analysis

  Swede.linmod = linmod(NextYear, LastYear, SwedeBetaList)

  Swede.ages = seq(0, 80, 2)
  Swede.beta1mat = eval.bifd(Swede.ages, Swede.ages, Swede.linmod$beta1estbifd)

# Figure 10.11

  persp(Swede.ages, Swede.ages, Swede.beta1mat,
        xlab="age", ylab="age",zlab="beta(s,t)",
        cex.lab=1.5,cex.axis=1.5)
}
# End Sweden example

##
## Section 10.5 Permutation Tests of Functional Hypotheses
##

#  Section 10.5.1 Functional t-Tests

#  ---------------  Comparing male and female growth data  ----------------

# Figure 10.12

ylim = with(growth, range(hgtm, hgtf))

with(growth, matplot(age, hgtm[, 1:10], type='l',
                     lty='dashed', ylab='height (cm)'))
with(growth, matlines(age, hgtf[, 1:10], lty='solid'))
legend('topleft', legend=c('girls', 'boys'),
       lty=c('solid', 'dashed'))

growthbasis = create.bspline.basis(breaks=growth$age, norder=6)
growfdPar = fdPar(growthbasis, 3, 10^(-0.5))

hgtffd = with(growth, smooth.basis(age,hgtf,growfdPar))
hgtmfd = with(growth, smooth.basis(age,hgtm,growfdPar))

tres = tperm.fd(hgtffd$fd,hgtmfd$fd)

# Figure 10.13

# Section 10.5.2 Functional F-Tests

#  ---------------  Testing for no effect of climate zone  ----------------

# temp36fd, regionList, betaList from Section 10.1.1 above

F.res = Fperm.fd(temp36fd, regionList, betaList)

# Figure 10.14
# plot in black and white

with(F.res,{
            q = 0.95
           ylims = c(min(c(Fvals, qval, qvals.pts)), max(c(Fobs,
                qval)))
            plot(argvals, Fvals, type = "l", ylim = ylims, col = 1,
                 lwd = 2, xlab = "day", ylab = "F-statistic",
                 cex.lab=1.5,cex.axis=1.5)
            lines(argvals, qvals.pts, lty = 3, col = 1, lwd = 2)
            abline(h = qval, lty = 2, col = 1, lwd = 2)
            legendstr = c("Observed Statistic", paste("pointwise",
                1 - q, "critical value"), paste("maximum", 1 -
                q, "critical value"))
            legend(argvals[1], 1.2, legend = legendstr,
                col = c(1, 1, 1), lty = c(1, 3, 2), lwd = c(2,
                  2, 2))
        }
)

##
## 10.6 Details for R Functions fRegress, fRegress.CV and fRegress.stderr
##
help(fRegress)
help(fRegress.CV)
help(fRegress.stderr)

##
## 10.7 Details for Function plotbeta
##
help(plotbeta)

##
## 10.8 Details for Function linmod
##
help(linmod)

##
## 10.9 Details for Functions Fperm.fd and tperm.fd
##
help(Fperm.fd)
help(tperm.fd)

##
## Section 10.10 Some Things to Try
##
# (exercises for the reader)

##
## Section 10.11  More to Read
##
