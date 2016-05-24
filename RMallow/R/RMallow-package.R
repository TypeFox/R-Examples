
#' Sample data set.
#' 
#' Simple synthetic data set containing 3 modal sequences in 15! space, with
#' some noise added.
#' 
#' @name datas
#' @docType data
#' @format The format is: num [1:1700, 1:15] 1 15 1 15 15 12 10 4 1 15 ...
#' @keywords datasets
#' @examples
#' data(datas)
#' head(datas)
NULL

#' 1980 APA Presidential Candidate ranking data.
#' 
#' This data is a pre-processed version of the 1980 American Psychological
#' Association Presidential candidate ranking data.  It has uninformative
#' rankings removed, and values pre-simplified into partial rankings.
#' 
#' 
#' @name elect
#' @docType data
#' @format The format is: int [1:1378, 1:3] 1 1 1 1 2 2 1 1 2 2 ...  - attr(*,
#' "dimnames")=List of 2 ..$ : chr [1:1378] "1" "2" "3" "6" ...  ..$ : chr
#' [1:3] "Carter" "Reagan" "Anderson"
#' @source The American Psychological Association, http://www.electionstudies.org/studypages/1980prepost/1980prepost.htm
#' @keywords datasets
#' @examples 
#' data(elect)
#' head(elect)
NULL





#' Fit Multi-modal Mallows' models to ranking data.
#' 
#' Fits the Mallows' model to ranking data.  Data can be partially or
#' fully-ranked.
#' 
#' \tabular{ll}{ Package: \tab RMallow\cr Type: \tab Package\cr Version: \tab
#' 1.0\cr Date: \tab 2012-02-18\cr License: \tab GPL (>= 2) } 
#' @name RMallow-package
#' @aliases RMallow-package RMallow
#' @docType package
#' @author Erik Gregory
#' Maintainer: <egregory2007@@yahoo.com>
#' @references "Mixtures of distance-based models for ranking data". Thomas 
#' Brendan Murphy & Donal Martin. 1 April 2002. Computational Statistics & 
#' Data Analysis 41 (2003) 645-655.
#'
#' @references "Estimating a Population Distribution of Sequences of k Items from Cross-
#' Sectional Data". Laurel A. Smith (Beckett) and Denis A. Evans. Journal of
#' the Royal Statistical Society.  Series C (Applied Statistics). Vol. 40,
#' No. 1 (1991), pp.31-42.  Blackwell Publishing for the Royal Statistical Society.
#' Accessed 16/08/2010. http://www.jstor.org/stable/2347903 .
#'
#' @references "A Non-iterative procedure for maximum likelihood estimation of the parameters of Mallows' Model Based on Partial Rankings".  Laura Adkins and Michael Flinger. Communication in Statistics - Theory and Methods, 27:9, 2199-2220. 1998, Marchel Dekker, Inc. http://dx.doi.org/10.1080/03610929808832223 .
#' @keywords ranking
NULL

#' Fitted version of the toy datas data set, with three modal sequences.
#' 
#' The data has 3 modal sequences, and we can compare this to the two.mode data
#' set.
#' 
#' 
#' @name three.mode
#' @docType data
#' @format The format is: List of 5 $ R :List of 3 ..$ : int [1:15] 1 2 3 4 5 6
#' 7 8 9 10 ...  ..$ : int [1:15] 1 3 5 7 9 2 4 6 8 10 ...  ..$ : int [1:15] 15
#' 14 13 12 11 10 9 8 7 6 ...  $ p : num [1:3] 0.447 0.118 0.435 $ lambda : num
#' [1:3] 2.01 1000 2.04 $ datas :'data.frame': 1700 obs. of 23 variables: ..$
#' X1 : num [1:1700] 1 15 1 15 15 12 10 4 1 15 ...  ..$ X2 : num [1:1700] 2 14
#' 2 14 14 13 13 12 2 14 ...  ..$ X3 : num [1:1700] 3 13 3 13 13 2 4 6 3 13 ...
#' ..$ X4 : num [1:1700] 4 12 4 12 12 8 7 1 4 12 ...  ..$ X5 : num [1:1700] 5
#' 11 5 11 11 9 14 5 5 11 ...  ..$ X6 : num [1:1700] 6 10 6 10 10 1 8 10 6 10
#' ...  ..$ X7 : num [1:1700] 7 9 7 9 9 15 1 13 7 9 ...  ..$ X8 : num [1:1700]
#' 8 8 8 8 8 10 9 9 8 8 ...  ..$ X9 : num [1:1700] 9 7 9 7 7 6 5 14 9 7 ...
#' ..$ X10 : num [1:1700] 10 6 10 6 6 11 11 8 10 6 ...  ..$ X11 : num [1:1700]
#' 11 5 11 5 5 3 15 2 11 5 ...  ..$ X12 : num [1:1700] 12 4 12 4 4 14 12 11 12
#' 4 ...  ..$ X13 : num [1:1700] 13 3 13 3 3 7 2 7 13 3 ...  ..$ X14 : num
#' [1:1700] 14 2 14 2 2 5 3 15 14 2 ...  ..$ X15 : num [1:1700] 15 1 15 1 1 4 6
#' 3 15 1 ...  ..$ clust : int [1:1700] 1 3 1 3 3 3 3 1 1 3 ...  ..$ pvals.1:
#' num [1:1700] 1.00 1.03e-91 1.00 2.04e-93 1.03e-91 ...  ..$ pvals.2: num
#' [1:1700] 0 0 0 0 0 0 0 0 0 0 ...  ..$ pvals.3: num [1:1700] 1.02e-92 1.00
#' 1.34e-93 1.00 1.00 ...  ..$ seq : Factor w/ 3 levels "1 2 3 4 5 6 7 8 9 10
#' 11 12 13 14 15",..: 1 3 1 3 3 3 3 1 1 3 ...  ..$ dists.1: num [1:1700] 0 105
#' 0 105 105 61 58 46 0 105 ...  ..$ dists.2: num [1:1700] 10 95 10 95 95 61 54
#' 54 10 95 ...  ..$ dists.3: num [1:1700] 105 0 105 0 0 44 47 59 105 0 ...  $
#' min.like: num [1:100] -122710 -51439 -50310 -49976 -49718 ...
#' @keywords datasets
#' @examples
#' data(three.mode)
#' head(three.mode[[4]])
NULL





#' Two-mode Mallows' model fit to toy data set "datas"
#' 
#' "datas" has 3 modes, but we observe here what happens when we try to fit it
#' with 2 modal sequences.  The most prominent modal sequences are 1:15, 15:1
#' 
#' 
#' @name two.mode
#' @docType data
#' @format The format is: List of 5 $ R :List of 2 ..$ : int [1:15] 1 2 3 4 5 6
#' 7 8 9 10 ...  ..$ : int [1:15] 15 14 13 12 11 10 9 8 7 6 ...  $ p : num
#' [1:2] 0.557 0.443 $ lambda : num [1:2] 2.05 2.02 $ datas :'data.frame': 1700
#' obs. of 21 variables: ..$ X1 : num [1:1700] 1 15 1 15 15 12 10 4 1 15 ...
#' ..$ X2 : num [1:1700] 2 14 2 14 14 13 13 12 2 14 ...  ..$ X3 : num [1:1700]
#' 3 13 3 13 13 2 4 6 3 13 ...  ..$ X4 : num [1:1700] 4 12 4 12 12 8 7 1 4 12
#' ...  ..$ X5 : num [1:1700] 5 11 5 11 11 9 14 5 5 11 ...  ..$ X6 : num
#' [1:1700] 6 10 6 10 10 1 8 10 6 10 ...  ..$ X7 : num [1:1700] 7 9 7 9 9 15 1
#' 13 7 9 ...  ..$ X8 : num [1:1700] 8 8 8 8 8 10 9 9 8 8 ...  ..$ X9 : num
#' [1:1700] 9 7 9 7 7 6 5 14 9 7 ...  ..$ X10 : num [1:1700] 10 6 10 6 6 11 11
#' 8 10 6 ...  ..$ X11 : num [1:1700] 11 5 11 5 5 3 15 2 11 5 ...  ..$ X12 :
#' num [1:1700] 12 4 12 4 4 14 12 11 12 4 ...  ..$ X13 : num [1:1700] 13 3 13 3
#' 3 7 2 7 13 3 ...  ..$ X14 : num [1:1700] 14 2 14 2 2 5 3 15 14 2 ...  ..$
#' X15 : num [1:1700] 15 1 15 1 1 4 6 3 15 1 ...  ..$ clust : int [1:1700] 1 2
#' 1 2 2 2 2 1 1 2 ...  ..$ pvals.1: num [1:1700] 1.00 4.15e-94 1.00 4.15e-94
#' 4.15e-94 ...  ..$ pvals.2: num [1:1700] 5.4e-93 1.0 5.4e-93 1.0 1.0 ...  ..$
#' seq : Factor w/ 2 levels "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15",..: 1 2 1 2 2
#' 2 2 1 1 2 ...  ..$ dists.1: num [1:1700] 0 105 0 105 105 61 58 46 0 105 ...
#' ..$ dists.2: num [1:1700] 105 0 105 0 0 44 47 59 105 0 ...  $ min.like: num
#' [1:100] -178063 -139298 -58290 -54074 -53902 ...
#' @keywords datasets
#' @examples
#' data(two.mode)
#' head(two.mode[[4]])
NULL





#' Bi-modal Mallow's model fit to the APA data set.
#' 
#' The two-modes seem to divide well between Democrats and Republicans...
#' 
#' @name two.seq
#' @docType data
#' @format The format is: List of 5 $ R :List of 2 ..$ : int [1:3] 1 3 2 ..$ :
#' int [1:3] 3 1 2 $ p : num [1:2] 0.541 0.459 $ lambda : num [1:2] 2.19 2.32 $
#' datas :'data.frame': 1378 obs. of 9 variables: ..$ Carter : int [1:1378] 1 1
#' 1 1 2 2 1 1 2 2 ...  ..$ Reagan : int [1:1378] 1 2 2 2 1 1 2 3 1 1 ...  ..$
#' Anderson: int [1:1378] 1 2 2 3 3 3 3 2 3 3 ...  ..$ clust : int [1:1378] 1 1
#' 1 1 2 2 1 1 2 2 ...  ..$ pvals.1 : num [1:1378] 0.541 0.992 0.992 0.932
#' 0.131 ...  ..$ pvals.2 : num [1:1378] 0.45893 0.00809 0.00809 0.06802
#' 0.86945 ...  ..$ seq : Factor w/ 2 levels "1 3 2","3 1 2": 1 1 1 1 2 2 1 1 2
#' 2 ...  ..$ dists.1 : num [1:1378] 0 0 0 1 2 2 1 0 2 2 ...  ..$ dists.2 : num
#' [1:1378] 0 2 2 2 1 1 2 3 1 1 ...  $ min.like: num [1:100] -6421 -3386 -2916
#' -2811 -2799 ...
#' @source American Psychological Association http://www.electionstudies.org/studypages/anes_mergedfile_1980/anes_mergedfile_1980.htm
#' @keywords dataset
#' @examples
#' data(two.seq)
#' head(two.seq[[4]])
NULL



