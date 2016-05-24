`inverseLogit` <-
function(x) {
#
# copyright 2006-present, University of Washington. All rights reserved.
# for terms of use, see the LICENSE file

              if (is.na(x)) return(NA)
              if (x >= 0) {
                if (-x >= log(.Machine$double.eps)) {
                  x <- exp(-x)
                  1/(1+x)
                }
                else 1
              }
             else {
                if (x >= log(.Machine$double.xmin)) {
                  x <- exp(x)
                  x/(1+x)
                }
                else 0
             }
            }

