#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# FUNCTION:                 DESCRIPTION:
#  rowCumsums,ANY            Computes cumulated sums by row
#  rowCumsums,timeSeries     Computes cumulated sums by row for timeSeries
################################################################################


setMethod("rowCumsums", "ANY", function(x, na.rm = FALSE, ...)
      {
          # Transform:
          if (!inherits(x, 'matrix'))
              x <- as(x, "matrix")

          if (na.rm)
              x <- na.omit(x)

          ans <- apply(x, 1, cumsum, ...)

          # special treatment when x has one row because apply returns a vector
          if (NCOL(x) > 1)
              t(ans)
          else
              matrix(ans, ncol = 1, dimnames = dimnames(x))
      })

      
# ------------------------------------------------------------------------------


setMethod("rowCumsums", "timeSeries", function(x, na.rm = FALSE, ...)
          setDataPart(x, callGeneric(getDataPart(x), na.rm = na.rm, ...)))

          
################################################################################

