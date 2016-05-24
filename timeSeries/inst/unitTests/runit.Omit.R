
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


test.naOmitMatrix =
function()
{
    x = as.timeSeries(data(LPP2005REC))[1:20, 1:4]
    colnames(x) = abbreviate(colnames(x), 6)

    x[1, 1] = NA
    x[3:4, 2] = NA
    x[18:20, 4] = NA
    show(x)

    timeSeries:::.naOmitMatrix(as.matrix(x))
    timeSeries:::.naOmitMatrix(as.matrix(x), "s")
    timeSeries:::.naOmitMatrix(as.matrix(x), "z")
    timeSeries:::.naOmitMatrix(as.matrix(x), "ir")
    timeSeries:::.naOmitMatrix(as.matrix(x), "iz")
    timeSeries:::.naOmitMatrix(as.matrix(x), "ie")

    # Return Value:
    return()

}


# ------------------------------------------------------------------------------


test.na.omit =
function()
{
    x = as.timeSeries(data(LPP2005REC))[1:20, 1:4]
    colnames(x) = abbreviate(colnames(x), 6)

    x[1, 1] = NA
    x[3:4, 2] = NA
    x[18:20, 4] = NA
    show(x)

    na.omit(x)
    na.omit(x, "s")
    na.omit(x, "z")
    na.omit(x, "ir")
    na.omit(x, "iz")
    na.omit(x, "ie")

    # Return Value:
    return()
}


################################################################################

