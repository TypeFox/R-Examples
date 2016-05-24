
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


test.signalSeries.Internal <-
function()
{
    # RUnit Test:

    x = rnorm(12)
    y = rnorm(12)
    timeSeries:::.signalSeries(as.matrix(x), 1:12)
    timeSeries:::.signalSeries(as.matrix(cbind(x,y)), 1:12)
}


# ------------------------------------------------------------------------------


test.timeSeries.Internal <-
function()
{
    # this is to test the problem when a ts object is passed to
    # timeSeries. It seems that as.matrix does not convert the object
    # to a matrix !!!
    z <- ts(matrix(rnorm(300), 100, 3), start=c(1961, 1), frequency=12)
    # class(as.matrix(z)) #<< mts ts and not matrix in R 2.9.0
    # Note that is is possible that a ts object is considered as a
    # matrix when timeSeries method as dispatched. Hence this check
    t <- timeSeries(z)
    checkTrue(identical(as(z, "matrix"), as(t, "matrix")))

}


################################################################################

