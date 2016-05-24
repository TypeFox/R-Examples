
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

# Copyrights (C)
# for this R-port:
#  1999 - 2007, Diethelm Wuertz, GPL
#  Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#  info@rmetrics.org
#  www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#  see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#  see Rmetrics's copyright file


################################################################################
# FUNCTION:                 DESCRIPTION:
#  rollFun                   Compute Rolling Function Value
#   rollMean                  Compute Rolling Mean
#   rollVar                   Compute Rolling Variance
#   rollMin                   Compute Rolling Minimum
#   rollMax                   Compute Rolling Maximum
################################################################################


test.rollingVector =
function()
{
    # Period:
    n = 3

    # TRIM = TRUE | na.rm = TRUE
    trim = TRUE
    na.rm = TRUE
    x = 1:10
    x[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = TRUE | na.rm = FALSE
    trim = TRUE
    na.rm = FALSE
    x = 1:10
    x[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    # ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    # print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = FALSE | na.rm = TRUE
    trim = FALSE
    na.rm = TRUE
    x = 1:10
    x[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = FALSE | na.rm = FALSE
    trim = FALSE
    na.rm = FALSE
    x = 1:10
    x[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    # ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    # print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.rollingTimeSeries =
function()
{
    # Time Series:
    charvec = paste("1999", 10, 11:20, sep = "-")
    print(charvec)
    ts = timeSeries(data = 1:10, charvec, units = "SERIES", zone = "GMT",
        FinCenter = "GMT")
    print(ts)

    # Period:
    n = 3

    # TRIM = TRUE | na.rm = TRUE
    trim = TRUE
    na.rm = TRUE
    x = ts
    series(x)[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = TRUE | na.rm = FALSE
    trim = TRUE
    na.rm = FALSE
    x = ts
    series(x)[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    # ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    # print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = FALSE | na.rm = TRUE
    trim = FALSE
    na.rm = TRUE
    x = ts
    series(x)[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # TRIM = FALSE | na.rm = FALSE
    trim = FALSE
    na.rm = FALSE
    x = ts
    series(x)[6] = NA
    cat("\ntrim: ", trim, "\n")
    cat("\n\nna.rm: ", na.rm, "\n")
    # Sum:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = sum)
    print(ans)
    # Mean:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = mean)
    print(ans)
    # Var:
    # ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = var)
    # print(ans)
    # Min:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = min)
    print(ans)
    # Max:
    ans = rollFun(x, n = n, trim = trim, na.rm = na.rm, FUN = max)
    print(ans)

    # Return Value:
    return()
}


################################################################################

