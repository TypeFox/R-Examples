
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


test.dim.timeSeries =
function()
{
    # RUnit Test:

    # Univariate Case:
    tS = timeSeries(format = "counts")[, 1]
    dim(tS)
    NCOL(tS)
    NROW(tS)
    ncol(tS)
    nrow(tS)
    dimnames(tS)
    colnames(tS)
    rownames(tS)
    # colnames<-
    # rownames<-
    is.array(tS)

    # Multivariate Case:
    tS = timeSeries(format = "counts")
    dim(tS)
    NCOL(tS)
    NROW(tS)
    ncol(tS)
    nrow(tS)
    dimnames(tS)
    colnames(tS)
    rownames(tS)
    # colnames<-
    # rownames<-
    is.array(tS)

    # Univariate Case:
    tS = timeSeries()[, 1]
    dim(tS)
    NCOL(tS)
    NROW(tS)
    ncol(tS)
    nrow(tS)
    dimnames(tS)
    colnames(tS)
    rownames(tS)
    # colnames<-
    # rownames<-
    is.array(tS)

    # Multivariate Case:
    tS = timeSeries()
    dim(tS)
    NCOL(tS)
    NROW(tS)
    ncol(tS)
    nrow(tS)
    dimnames(tS)
    colnames(tS)
    rownames(tS)
    # colnames<-
    # rownames<-
    is.array(tS)
}


################################################################################

