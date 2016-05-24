
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


test.signalCounts <-
function()
{
    # RUnit Test:

    int = c(1, 10, 100, 21, 135)
    print(timeSeries:::.signalCounts(sample(int)))

    nc = timeSeries:::.signalCounts(int)
    nc

    ns = sample(nc)
    ns

    sorted = sort(ns)
    sorted
    as.integer(sorted)
    ns

    ordered = order(ns)
    ordered
    ns[ordered]
    as.integer(ns[ordered])

    timeSeries:::.signalCounts(1:12)
    timeSeries:::.signalCounts(sample(1:12))
    timeSeries:::.signalCounts(timeSeries:::.signalCounts(1:12))
}


################################################################################

