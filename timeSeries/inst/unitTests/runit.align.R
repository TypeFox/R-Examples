
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


test.align.timeSeries <-
function()
{
    # RUnit Test:

    # .align.timeSeries(x, method = c("before", "after", "interp"),
    #   startOn = "hours", by = "30 m")

    set.seed(1953)
    tD = timeCalendar(
        y = rep(2008, times = 6), m = rep(4, times = 6), d = rep(10:11, each = 3),
        h = sample(1:23)[1:6], min = sample(1:59)[1:6], s = sample(1:59)[1:6])
    tS = timeSeries(rnorm(6), tD)
    align(tS)
    align(tS, method="interp")

    # Note, we should als add an argument to trim NAs
}

################################################################################

