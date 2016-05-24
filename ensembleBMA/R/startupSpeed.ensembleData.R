startupSpeed.ensembleData <-
function (x) 
{
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file
#
 if (is.null(startup <- attr( x, "startupSpeed"))) {
   class(x) <- "data.frame" 
   x$startup
 }
 else startup
}

