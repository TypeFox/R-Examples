#  -----------------------------------------------------------------------
#                            Growth Data Analyses
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#
#                          Overview of the analyses
#
#  These analyses are intended to illustrate the analysis of nonperiod data
#  where a spline basis is the logical choice.  These analyses complement
#  the daily growth data in that sense.
#  The growth data have the additional feature of being essentially
#  monotonic or, to say the same thing in another way, have an essentially
#  positive first derivative or velocity.  This requires monotone smoothing.
#  Moreover, most of the interpretability of the growth data comes from
#  inspecting the acceleration of the height curves, so that great emphasis
#  is placed here on getting a good sensible and stable acceleration
#  estimate.
#  Finally, a large prortion of the variation in the growth curve data is due
#  to phase variation, mainly through the variation in the timing of the
#  pubertal growth spurt.  Registration therefore plays a major role and is
#  especially illustrated here.
#  Most of the analyses are carried out on the Berkeley growth data, which
#  have the advantage of being freely distributable, whereas as more recent
#  and larger data bases require special permission from the agencies that
#  are responsible for them.  Not much is lost, however, since the quality
#  of the Berkeley data are quite comparable to those of other datasets.
#  The primary analyses are the monotone smoothing of the data.  The right
#  smoothing level is taken as known, and was determined by other analyses
#  in the Matlab language.  The monotone smoothing function used here
#  requires the use of low-level code in C and C++, but even with that help,
#  computation times are substantially longer than in Matlab.
#  Following monotone smoothing, the growth data are registered, an
#  essential step because of the large variation in the timing of the
#  pubertal growth spurt.  The pubertal growth spurts are aligned using
#  landmark registration, and the land-mark registered curves are then
#  registered using continuous registration.
#  The final analysis is of a set of data on a single boy where the
#  measurements are taken every three days or so, rather than twice a year.
#  These data show that growth is rather more complex than the traditional
#  data could have revealed.
#  -----------------------------------------------------------------------

#  -----------------------------------------------------------------------
#                           Berkeley Growth Data
#  -----------------------------------------------------------------------

#  Last modified 9 March 200

#  ------------------------  input the data  -----------------------

age  <- c( seq(1, 2, 0.25), seq(3, 8, 1), seq(8.5, 18, 0.5))
nage <- length(age)

ncasem <- 39
ncasef <- 54

hgtm <- t(matrix(scan("../data/hgtm.txt",0),ncasem,nage,byrow=T))
hgtf <- t(matrix(scan("../data/hgtf.txt",0),ncasef,nage,byrow=T))

#  Save the data

growthdata <- list(hgtm  = hgtm,  hgtf = hgtf, age = age)

save(growthdata, file = "growthdata")

