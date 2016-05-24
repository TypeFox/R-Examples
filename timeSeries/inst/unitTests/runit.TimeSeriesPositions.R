
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


test.time =
function()
{
    # Generate nivariate daily random sequence
    set.seed(4711)
    data = round(rnorm(12), 2)
    charvec = timeCalendar(2006)
    uTS = timeSeries(data, charvec, units = "uTS")
    uTS

    # Get Positions:
    POS = time(uTS)
    POS
    checkIdentical(charvec, POS)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


"test.time<-" =
function()
{
    # Generate nivariate daily random sequence
    set.seed(4711)
    data = round(rnorm(12), 2)
    charvec = timeCalendar(2006)
    uTS = timeSeries(data, charvec, units = "uTS")
    uTS

    # Add one Day to Positions:
    POS = time(uTS)
    time(uTS) <- POS + 24*3600
    uTS

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.timeSeriesOrdering =
function()
{
    #  sample.timeSeries - Resamples a 'timeSeries' object in time
    #  sort.timeSeries - Sorts reverts a 'timeSeries' object in time
    #  rev.timeSeries - Reverts a 'timeSeries' object in time
    #  start.timeSeries - Extracts start date of a 'timeSeries' object
    #  end.timeSeries - Extracts end date of a 'timeSeries' object

    # Generate univariate monthly random sequence:
    set.seed(4711)
    data = cbind(1:12, round(rnorm(12), 2))
    positions = timeCalendar(2006)
    uTS = timeSeries(data, positions)
    uTS

    # Sample/Sort:
    target = uTS
    target
    # current = sort(sample(uTS))
    # current
    # checkIdentical(target, current)

    # Revert:
    target = uTS
    target
    current = rev(rev(uTS))
    current
    checkTrue(!sum(target - current))

    # Start/End date of Series:
    start(uTS)
    end(uTS)

    # Return Value:
    return()
}


################################################################################

