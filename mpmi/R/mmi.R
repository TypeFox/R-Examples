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


mmi <-
function(cts, disc, level = 3L, na.rm = FALSE, h, ...)
{
    cts <- as.matrix(cts)
    mode(cts) <- "double"

    mi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))
    bcmi <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))
    zans <- matrix(0.0, nrow = ncol(cts), ncol = ncol(disc))

    # Convert to integers (probably slow)
    dint <- matrix(0L, nrow = nrow(disc), ncol = ncol(disc))
    for (i in 1: ncol(disc))
    {
        dint[,i] <- as.integer(factor(disc[,i]))
    }

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
 
    result <- .Fortran("mmim", 
                    cts, 
                    as.integer(nrow(cts)), 
                    as.integer(ncol(cts)),
                    dint,
                    as.integer(nrow(dint)), 
                    as.integer(ncol(dint)),
                    mi = mi, 
                    bcmi = bcmi, 
                    zvalues = zans,
                    as.double(h), 
                    NAOK = TRUE, 
                    DUP = TRUE)

    return(result[c("mi", "bcmi", "zvalues")])
}

mmi.pw <- function(cts, disc, h, ...)
{
    if (length(cts) != length(disc)) stop("Input vectors must be the same length")

    # Remove missing values pairwise
    ok <- !is.na(disc) & !is.na(cts)
    disc <- disc[ok]
    cts <- cts[ok]

    if (missing(h))
    {
        h <- dpikSafe(cts, level = 3L, kernel = "epanech", ...)
    }

    lf <- length(cts)

    result <- .Fortran("mmipw", cts = as.double(cts), 
                        lc = as.integer(lf), 
                        disc = as.integer(factor(disc)),
                        h = as.double(h),
                        mi = as.double(0), 
                        bcmi = as.double(0),
                        zvalue = as.double(0),
                        DUP = TRUE)

    return(result[c("mi", "bcmi", "zvalue")])
}

