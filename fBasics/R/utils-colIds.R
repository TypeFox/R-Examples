
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
# GENERATION:               DESCRIPTION:
#  colIds                    Retrieves column names of a matrix-like object
#  rowIds                    Retrieves row names of a matrix-like object
#  colIds<-                  Sets column names of a matrix-like object
#  rowIds<-                  Sets row names of a matrix-like object
################################################################################


colIds <-
function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Retrieves row names of a matrix-like object

    # FUNCTION:

    # Convert to Matrix
    x = as.matrix(x)

    # Return Value:
    colnames(x, ...)
}


# ------------------------------------------------------------------------------


rowIds <-
function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Retrieves row names of a matrix-like object

    # FUNCTION:

    # Convert to Matrix
    x = as.matrix(x)

    # Return Value:
    rownames(x, ...) }


# ------------------------------------------------------------------------------


"colIds<-" <-
function(x, value)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets column names of a matrix-like object

    # FUNCTION:

    # Column Names:
    dn = dimnames(x)
    if(is.null(dn)) {
        if(is.null(value)) return(x)
        if((nd = length(dim(x))) < 2)
            stop("Object has less than two dimensions")
        dn = vector("list", nd)
    }
    if(length(dn) < 2)
        stop("Object has less than two dimensions")
    if(is.null(value)) dn[2] = list(NULL) else dn[[2]] = value
    dimnames(x) = dn

    # Return Value:
    x
}


# ------------------------------------------------------------------------------


"rowIds<-" <-
function(x, value)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Sets row names of a matrix-like object

    # FUNCTION:

    # Row names:
    dn = dimnames(x)
    if(is.null(dn)) {
        if(is.null(value)) return(x)
        if((nd = length(dim(x))) < 1)
            stop("attempt to set rownames on object with no dimensions")
        dn = vector("list", nd) }
    if(length(dn) < 1)
        stop("attempt to set rownames on object with no dimensions")
    if(is.null(value)) dn[1] = list(NULL) else dn[[1]] = value
    dimnames(x) = dn

    # Return Value:
    x
}


################################################################################

