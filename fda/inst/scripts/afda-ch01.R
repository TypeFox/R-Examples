###
###
### Ramsey & Silverman (2002) Applied Functional Data Analysis (Springer)
###
### ch. 1.  Introduction 
###
library(fda)
##
## Intro to ch. 2.  Criminology
##
# p. 3, Figure 1.1.  An individual in the criminology sample

# Data not available.  







# p. 3-4, Figure 1.2.  413 subjects in the criminology study

# Data not available.  







##
## Intro to ch. 3.  Nondurable goods index 
##
# pp.  4-6, Figure 1.3.  US nondurable goods index 1919-2000
#durtime = (0:(ndur-1))./12 + 1919;
ndur <- length(nondurables)
durtime = (0:(ndur-1))/12 + 1919;
lognondur = log10(nondurables);

plot(nondurables, xlab="Year", ylab="Nondurable goods index")
plot(nondurables, log="y", xlab="Year",
     ylab="Nondurable goods index")

# pp. 5-6, Figure 1.4.  Phase-plane plots for 1923 & 1996
#  smooth the log data with order 8 splines, knots at data points

#  Fit smooth per sec. 3.6.
goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)

LfdobjNonDur     = int2Lfd(4);

#goodsfdPar = fdPar(goodsbasis, LfdobjNonDur, lambda=1e-11)
#lognondursmth = smooth.basis(durtime, coredata(lognondur), goodsfdPar);
logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)

#str(lognondursmth)

nondur1964.1967 <- window(nondurables, 1964, 1967)

plot(log10(nondur1964.1967), type="p", axes=FALSE, xlab="Year",
     ylab=expression(paste(log[10], " nondurable goods index")) )
axis(2)
axis(1, 1964:1967)
axis(1, seq(1964, 1967, by=0.5), labels=FALSE)

#durtimefine = linspace(1964,1967,101);
durtimefine <- seq(1964, 1967, length=181)

#fit = eval.fd(durtimefine, lognondursmth);
logNondurSm1964.67 = eval.fd(durtimefine, logNondurSm$fd);
lines(durtimefine, logNondurSm1964.67)
abline(v=1965:1966, lty=2)

##
## Intro to ch. 4.  Bone shapes in paleopathology 
##
# pp. 6-7, Figure 1.5.  Digital image of a femur
# Data not available.  

                                         









##
## Intro to ch. 5.  ADHD Reaction time
##
# pp. 7-8, Figure 1.6.  Reaction time distributions for two children 
# Data not available.











##
## Intro to ch. 6.  Human growth
##
# pp.  8-9, Figure 1.7.  Raw growth data for one individual 
with(growth, plot(age, hgtm[, 1], pch="+",
                  ylab="Measured height (cm.)"))

##
## Intro to ch. 7.  Time warping handwriting and weather
##
# pp. 9-10, Figure 1.8, "fda" written by hand 20 times
# Data not available.

# 'handwrit' is cursive;  Figure 1.8 is print.  
plot(handwrit[, 1,], type="l")

#?????????????????????









##
## Intro to ch. 8.  Bone shapes and arthritis
##
# no data used in this intro

##
## Intro to ch. 9.  Test Items
##
# pp. 11-12.  Figure 1.9.  Probability of success on two test questions
# Data not available.  










##
## Intro to ch. 10.  Lip acceleration
##
# no data used in this intro

##
## Intro to ch. 11.  Handwriting printed characters
##
# no data used in this intro

##
## Intro to ch. 12.  Juggling
##
# p. 14.  Figure 1.10.  A juggling cycle 
# Data not available.  










