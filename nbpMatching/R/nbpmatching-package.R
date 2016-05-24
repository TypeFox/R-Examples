#'Nonbipartite Matching
#'
#'This package will take an input distance matrix and generate the set of
#'pairwise matches that minimizes the sum of distances between the pairs by
#'running nonbimatch.
#'
#'The most current documentation is available at
#'\url{http://biostat.mc.vanderbilt.edu/wiki/Main/MatchedRandomization}.
#'
#'
#'@name nbpMatching-package
#'@aliases nbpMatching-package nbpMatching
#'@docType package
#'@author Bo Lu, Robert Greevy, Cole Beck
#'
#'Maintainer: Cole Beck \email{cole.beck@@vanderbilt.edu}
#'@references Lu B, Greevy R, Xu X, Beck C. Optimal Nonbipartite Matching and
#'its Statistical Applications. The American Statistician. Vol. 65, no. 1. :
#'21-30. 2011.
#'
#'Greevy RA Jr, Grijalva CG, Roumie CL, Beck C, Hung AM, Murff HJ, Liu X,
#'Griffin MR. Reweighted Mahalanobis distance matching for cluster-randomized
#'trials with missing data. Pharmacoepidemiol Drug Saf. 2012 May;21 Suppl
#'2:148-54. doi: 10.1002/pds.3260.
#'@keywords package cluster array
#'@import methods
#'@import stats
#'@importFrom Hmisc transcan hdquantile
#'@importFrom MASS ginv
#'@importFrom utils read.csv
#'@useDynLib nbpMatching mwrap
#'@examples
#'
#'# create a covariate matrix
#'df <- data.frame(id=LETTERS[1:25], val1=rnorm(25), val2=rnorm(25))
#'# create distances
#'df.dist <- gendistance(df, idcol=1)
#'# create distancematrix object
#'df.mdm <- distancematrix(df.dist)
#'# create matches
#'df.match <- nonbimatch(df.mdm)
#'# review quality of matches
#'df.qom <- qom(df.dist$cov, df.match$matches)
#'
#'# some helper functions are available
#'# runner -- start with the covariate, run through the entire process
#'df.1 <- runner(df, idcol=1)
#'# full.qom -- start with the covariate, generate a full quality of match report
#'df.2 <- full.qom(df)
#'
#'\dontrun{
#'try a large matrix
#'nonbimatch(distancematrix(as.matrix(dist(sample(1:10^8, 5000, replace=TRUE)))))
#'}
#'
NULL

#'Internal nbpMatching objects.
#'
#'Internal nbpMatching objects.
#'
#'This function should not be called by the user.
#'
#'@name nbpMatching-internal
#'@aliases .requireCachedGenerics initialize,distancematrix-method
#'[<-,distancematrix,ANY,ANY,ANY-method [<-,distancematrix-method
#'[[<-,distancematrix,ANY,ANY,ANY-method [[<-,distancematrix-method
#'runner runner,data.frame-method full.qom full.qom,data.frame-method
#'show,nonbimatch-method
#'@exportMethod runner
#'@exportMethod full.qom
#'@keywords internal
#'
NULL

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("Notice:
Formerly the gendistance() function scaled the Mahalanobis distances into large
integers, as required by the nonbimatch() function. Starting in version 1.5.0,
gendistance() will return unscaled distances. This facilitates comparison to an
appropriate F distribution for multivariate normal data. Any required scaling
will happen invisibly within nonbimatch(). This notice will be removed in a
future version of nbpMatching.")
}
