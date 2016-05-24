`ensembleForecasts.ensembleData` <-
function (x) 
{ 
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 k <- attr(x, "ensembleSize")
 attr(x, "ensembleSize") <- NULL
 class(x) <- "data.frame"
 as.matrix(x[, 1:k])
}

