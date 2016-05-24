
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


test.as <-
function()
{
    # RUnit Test:

    # Note, you can also use ...
    is(timeSeries(), "timeSeries")

    # Series
    ts = dummySeries()
    x = timeSeries:::.signalSeries(as.matrix(ts))
    y = timeSeries:::.timeSeries(as.matrix(ts), as.numeric(time(ts), "sec"))

    # A vector to a timeSeries
    as.vector(x)
    as.vector(x[,1])
    as.vector(y)
    as.vector(y[,1])

    # as.numeric:
    as.numeric(x)
    as.numeric(x[,1])
    as.numeric(y)
    as.numeric(y[,1])
}


################################################################################

