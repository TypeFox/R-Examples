
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


test.removeNA =
function()
{

    # Create matrix object:
    set.seed(1985)
    M = 5
    N = 20
    x = matrix(round(rnorm(M*N), 3), ncol = M)
    colnames(x) = 1:M
    rownames(x) = 1:N
    nNA = 10
    nCol = trunc(runif(nNA, 1, M+1))
    nRow = trunc(runif(nNA, 1, N+1))
    for (i in 1:nNA) x[nRow[i], nCol[i]] = NA
    print(x)
    ans = removeNA(x)
    print(ans)

    # Create data.frame object:
    x.df = as.data.frame(x)
    class(x.df)
    ans = removeNA(x.df)
    print(ans)
    class(ans)

    # Create timeSeries object:
    tD = timeCalendar(m = 1, d = 1:N)
    x.tS = timeSeries(x, tD)
    print(x.tS)
    ans = removeNA(x.tS)
    print(ans)
    class(ans)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.substituteNA =
function()
{
    # Create matrix object:
    set.seed(1985)
    M = 5
    N = 20
    x = matrix(round(rnorm(M*N), 3), ncol = M)
    colnames(x) = 1:M
    rownames(x) = 1:N
    nNA = 10
    nCol = trunc(runif(nNA, 1, M+1))
    nRow = trunc(runif(nNA, 1, N+1))
    for (i in 1:nNA) x[nRow[i], nCol[i]] = NA
    print(x)

    # Substitute:
    ans = substituteNA(x)
    print(ans)
    ans = substituteNA(x, "mean")
    print(ans)
    ans = substituteNA(x, "median")
    print(ans)

    # Create data.frame object:
    x.df = as.data.frame(x)
    print(x.df)
    class(x.df)

    # Substitute:
    ans = substituteNA(x.df)
    print(ans)
    ans = substituteNA(x.df, "mean")
    print(ans)
    ans = substituteNA(x.df, "median")
    print(ans)

    # Create timeSeries object:
    tD = timeCalendar(m = 1, d = 1:N)
    x.tS = timeSeries(x, tD)
    print(x.tS)
    class(x.tS)

    # Substitute:
    ans = substituteNA(x.tS)
    print(ans)
    ans = substituteNA(x.tS, "mean")
    print(ans)
    ans = substituteNA(x.tS, "median")
    print(ans)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.interpNA =
function()
{
    # Interpolate Column-by-Column

    # Create matrix object:
    set.seed(1985)
    M = 5
    N = 20
    x = matrix(round(rnorm(M*N), 3), ncol = M)
    colnames(x) = 1:M
    rownames(x) = 1:N
    nNA = 10
    nCol = trunc(runif(nNA, 1, M+1))
    nRow = trunc(runif(nNA, 1, N+1))
    for (i in 1:nNA) x[nRow[i], nCol[i]] = NA
    print(x)

    # Interpolate:
    ans = interpNA(x, "linear")
    print(ans)
    ans = interpNA(x, "before")
    print(ans)
    ans = interpNA(x, "after")
    print(ans)

    # Return Value:
    return()
}


################################################################################

