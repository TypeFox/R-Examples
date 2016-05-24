
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                      DESCRIPTION:
#  isDaily,timeSeries-method      Tests if a time series is a daily series
#  isMonthly,timeSeries-method    Tests if a time series is a monthly series
#  isQuarterly,timeSeries-method  Tests if a time series is a quarterly series
#  isRegular,timeSeries-method    Tests if a time series is a regular series
#  frequency,timeSeries-method    Returns the frequency of a regular time series
################################################################################


setMethod("isDaily", "timeSeries", function(x) callGeneric(time(x)))


setMethod("isQuarterly", "timeSeries", function(x) callGeneric(time(x)))


setMethod("isMonthly", "timeSeries", function(x) callGeneric(time(x)))


setMethod("isRegular", "timeSeries", function(x) callGeneric(time(x)))


setMethod("frequency", "timeSeries", function(x, ...) callGeneric(time(x), ...))


################################################################################

