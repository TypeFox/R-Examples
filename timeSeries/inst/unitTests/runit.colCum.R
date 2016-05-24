
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


test.colCum <-
    function()
{
    # RUnit Test:

    # Signal Series
    ts <- dummySeries(format = "counts")
    colCumsums(ts)
    colCummaxs(ts)
    colCummins(ts)
    colCumprods(ts)
    colCumreturns(ts)

    # Time Series:
    ts <- dummySeries()
    colCumsums(ts)
    colCummaxs(ts)
    colCummins(ts)
    colCumprods(ts)
    colCumreturns(ts)

    # check that timeSeries with one row still works ...
    t <- ts[1,]

    checkTrue(is(colCumsums(t), "timeSeries"))
    checkTrue(is(colCummaxs(t), "timeSeries"))
    checkTrue(is(colCummins(t), "timeSeries"))
    checkTrue(is(colCumprods(t), "timeSeries"))
    checkTrue(is(colCumreturns(t), "timeSeries"))

    checkEquals(nrow(colCumsums(t)), 1)
    checkEquals(nrow(colCummaxs(t)), 1)
    checkEquals(nrow(colCummins(t)), 1)
    checkEquals(nrow(colCumprods(t)), 1)
    checkEquals(nrow(colCumreturns(t)), 1)

}


################################################################################

