
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

test.na.contiguous =
    function()
{
    ## Dummy timeSeries with NAs entries
    data1 <- matrix(c(NA, 1), ncol = 2)
    data2 <- matrix(rep(2, 4), ncol = 2)
    data3 <- matrix(c(NA, 3), ncol = 2)
    data4 <- matrix(rep(4, 4), ncol = 2)

    data <- rbind(data1, data2, data3, data4)

    ts <- timeSeries(data, timeCalendar()[1:6])
    ## Find the longest consecutive non-missing values
    ans <- na.contiguous(ts)
    check <- getDataPart(ans)
    dimnames(check) <- NULL

    checkIdentical(data2, getDataPart(check))
}

