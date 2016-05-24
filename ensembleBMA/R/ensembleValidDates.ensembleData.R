ensembleValidDates.ensembleData <-
function (x) 
{ 
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 class(x) <- "data.frame"
 as.character(x$dates)
}

