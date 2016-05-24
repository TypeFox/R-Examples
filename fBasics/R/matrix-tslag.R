
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# TIME SERIES               DESCRIPTION:
#  tslag                     Lagged/leading vector/matrix of selected orders
#  .tslag1                   Internal Function used by tslag
################################################################################


tslag <-
function(x, k = 1, trim = FALSE)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a lagged or leading vector/matrix of selected order(s).

    # Arguments:
    #   x - a vector of data, missing values (NA) are allowed.
    #   k - the number of positions the new series is to lag
    #       or to lead the input series.
    #   trim - a logical flag, if TRUE, the missing values at the
    #       beginning or end of the returned series will be trimmed.
    #       The default value is FALSE.

    # Details:
    #   With a positive value of "k" we get a lagged series and with
    #   a negative value we get a leading series.

    # Examples:
    #   tslag(rnorm(10), 2)
    #   tslag(rnorm(10), -2:2)
    #   tslag(rnorm(10), -2:2, trim = TRUE)

    # FUNCTION:

    # Bind:
    ans = NULL
    for ( i in k) {
        ans = cbind(ans, .tslag1(x, i))
    }

    # Trim:
    if (trim) {
        indexes = (1:length(ans[,1]))[!is.na(apply(ans, 1, sum))]
        ans = ans[indexes, ]
    }

    # As Vector:
    if (length(k) == 1) ans = as.vector(ans)

    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.tslag1 <-
function(x, k)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Internal Function used by function tslag.

    # FUNCTION:
    y = x
    if (k > 0) y = c(rep(NA, times = k), x[1:(length(x)-k)])
    if (k < 0) y = c(x[(-k+1):length(x)], rep(NA, times = -k))

    # Return Value:
    y
}


################################################################################

