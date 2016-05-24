
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


test.lag <-
function()
{
    # RUnit Test:

    tS = round(dummySeries(flormat = "counts"), 3)[, 1]
    tS
    lag(tS)
    lag(tS, k = -2:2)
    lag(tS, k = -2:2, trim = TRUE)

    tS = round(dummySeries(), 3)[, 1]
    tS
    lag(tS)
    lag(tS, k = -2:2)
    lag(tS, k = -2:2, trim = TRUE)


    # check colnames when using multiple lag indexes.
    data <- matrix(runif(12), ncol = 2)
    charvec <- rev(paste("2009-0", 1:6, "-01", sep = ""))
    S <- timeSeries(data, charvec)
    colnames(S) <- paste("S", 1:2, sep = ".")
    ts <- lag(S, -1:1)
    checkIdentical(colnames(ts), c("S.1[-1]", "S.1[0]", "S.1[1]", "S.2[-1]",
                                   "S.2[0]", "S.2[1]"))

}



################################################################################

