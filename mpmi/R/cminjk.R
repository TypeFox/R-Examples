#  
#     Copyright 2013 Chris Pardy <cpardy@unsw.edu.au>
# 
#     This file is part of the mpmi R package.
# 
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, version 3.
# 
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
# 
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <http://www.gnu.org/licenses/>.
#  


cminjk <-
function(cts, level = 3L, na.rm = FALSE, h, ...)
{
    cts <- as.matrix(cts)
    mode(cts) <- "double"

    mi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(cts))
    # bcmi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(cts))
    # zans <- matrix(0.0, nrow = ncol(cts), ncol = ncol(cts))

    # Calculate bandwidths
    if (missing(h))
    {
        if(na.rm)
        {
            h2 <- function(x)
            {
                return(dpikSafe(x[is.finite(x)], level = level, kernel = "epanech", ...))
            }
            h <- apply(cts, 2, h2)
        } else 
        {
            h <- apply(cts, 2, dpikSafe, level = level, kernel = "epanech", ...)
        }
    }
 
    result <- .Fortran("cmimnjk", 
                    cts, 
                    as.integer(nrow(cts)), 
                    as.integer(ncol(cts)),
                    mi = mi, 
                    # bcmi = bcmi, 
                    # zvalues = zans,
                    as.double(h), 
                    NAOK = TRUE, 
                    DUP = TRUE)

    return(result$mi)
}

# Pairwise only
cminjk.pw <- function(v1, v2, h, ...)
{
    if (length(v1) != length(v2)) stop("Vectors must be of the same length")
    ok1 <- !is.na(v1)
    ok2 <- !is.na(v2)
    
    # Remove samples with missing values
    ok <- ok1 & ok2
    v1 <- v1[ok]
    v2 <- v2[ok]
    if (missing(h))

    {
        h <- c(dpikSafe(v1, level = 3L, kernel = "epanech", ...), dpikSafe(v2, level = 3L, kernel = "epanech", ...))
    }

    result <- .Fortran("cmipwnjk",
            as.double(v1),
            as.double(v2),
            as.integer(length(v1)),
            as.double(h[1]),
            as.double(h[2]),
            mi = as.double(0), 
            # bcmi = as.double(0),
            # zvalue = as.double(0),
            DUP = TRUE)

    return(result$mi)
}

