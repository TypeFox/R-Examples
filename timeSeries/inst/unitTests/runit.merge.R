
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


test.merge.timeSeries =
function()
{
    # RUnit Test:

    # Time Stamps:
    x = timeSeries()[,1]
    x
    y = timeSeries()
    y
    merge(x, y)

    # Signal Counts:
    x = timeSeries(format = "counts")[,1]
    x
    y = timeSeries(format = "counts")
    y
    merge(x, y)


    x <- dummySeries()[,1]
    x
    y <- dummySeries()
    y
    merge(x, y)


    # check that merge method can deal with timeSeries that have
    # colnames that are invalid data.frame colnames. For example
    # "S[-1]".

    data <- matrix(runif(18), ncol = 3)
    charvec <- rev(paste("2009-0", 1:6, "-01", sep = ""))
    S <- timeSeries(data, charvec)
    colnames(S) <- paste("S", 1:3, sep = ".")
    ts <- merge(S[,2], lag(S[,1], -1:1))

    checkIdentical(dim(ts), c(6L,4L))

}



################################################################################

