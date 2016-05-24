
# Rmetrics is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# Rmetrics is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################


test.subset <-
function()
{


    ts <- dummySeries()
    mat <- as.matrix(ts)

    # we want the same subset-ting rules as for a matrix
    # but we always print result in vertical style !

    # --------------------------------------------------------------------------
    # index

    checkIdentical(
                   ts[],
                   ts)

    checkTrue(suppressWarnings(is.na(ts[""])))
    checkTrue(is.na(mat[""]))

    checkIdentical(
                   as.matrix(ts[seq(4),2]),
                   mat[seq(4),2,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[rep(FALSE, 3), 1]),
                   mat[rep(FALSE, 3), 1,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[FALSE, 1]),
                   mat[FALSE, 1, drop = FALSE])

    checkIdentical(
                   as.matrix(ts[rep(TRUE), 2]),
                   mat[rep(TRUE), 2, drop=FALSE])

    charvec <- as.character(timeCalendar()[1:3])
    checkIdentical(
                   as.matrix(ts[charvec, 1]),
                   mat[charvec, 1, drop = FALSE])

    checkIdentical(
                   as.matrix(ts[seq(4),]),
                   mat[seq(4),,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[rep(FALSE, 3), ]),
                   mat[rep(FALSE, 3), ,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[FALSE, ]),
                   mat[FALSE, ,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[rep(TRUE), ]),
                   mat[rep(TRUE), ,drop=FALSE ])

    dd <- as.character(time(ts)[1])
    checkIdentical(
                   as.matrix(ts[dd, ]),
                   mat[dd, ,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[,2]),
                   mat[,2,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[2,FALSE]),
                   mat[2,FALSE, drop=FALSE])

    # prefer to have an empty timeSeries instead of empty data with row names
    checkIdentical(
                   as.matrix(ts[,FALSE]),
                   mat[,FALSE, drop = FALSE])

    checkIdentical(
                   as.matrix(ts[,TRUE ]),
                   mat[,TRUE ,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[, "TS.1"]),
                   mat[, "TS.1", drop = FALSE])


    # --------------------------------------------------------------------------
    # timeDate

    checkIdentical(
                   ts[timeCalendar()[1:5], 2],
                   ts[1:5,2])

    checkIdentical(
                   ts[timeCalendar()[1:5], ],
                   ts[1:5,])

    # --------------------------------------------------------------------------
    # logical matrix and timeSeries

    i <- ts < 0.4

    checkException(ts[series(i), ], silent = TRUE)
    checkException(ts[i, ], silent = TRUE)
    checkException(mat[series(i), ], silent = TRUE) # it fails as expected

    checkIdentical(
                   as.matrix(ts[series(i)[,1], ]),
                   mat[series(i)[,1], , drop=FALSE])

    checkIdentical(
                   as.matrix(ts[i[,1], ]),
                   mat[series(i)[,1], , drop=FALSE])

    checkIdentical(
                   as.matrix(ts[series(i)[,1],1]),
                   mat[series(i)[,1],1,drop=FALSE])

    checkIdentical(
                   as.matrix(ts[i[,1],1]),
                   mat[series(i)[,1],1,drop=FALSE])

    # this should fail
    checkException(ts[series(i), 2], silent = TRUE)
    checkException(ts[i, 2], silent = TRUE)
    checkException(ts[series(i), 1], silent = TRUE)

    checkException(ts[series(i),1], silent = TRUE)
    checkException(ts[i,1], silent = TRUE)
    checkException(mat[series(i),1], silent = TRUE)

    checkException(ts[series(i),], silent = TRUE)
    checkException(mat[series(i),], silent = TRUE)

    checkIdentical(
                    ts[series(i)],
                    mat[series(i)])

    checkIdentical(
                   ts[i],
                   mat[series(i)])


    # --------------------------------------------------------------------------
    # $,timeSeries method
    df <- as.data.frame(ts)

    checkIdentical(
                   ts$TS.,
                   df$TS.)

    checkIdentical(
                   ts$TS.1,
                   df$TS.1)

    checkIdentical(
                   ts$a,
                   df$a)

    colnames(ts) <- c("aa", "bb")
    colnames(df) <- c("aa", "bb")

    checkIdentical(
                   ts$a,
                   df$a)

    checkIdentical(
                   ts$b,
                   df$b)

}


################################################################################

