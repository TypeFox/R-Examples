
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


test.order <-
function()
{
    # RUnit Test:

    tS = timeSeries(matrix(rnorm(500), 100), units = sample(LETTERS[1:5]))
    head(tS)
    colnames(tS)
    sortColnames(tS)
    sampleColnames(tS)
    orderColnames(tS)
    statsColnames(tS, FUN = colMeans)
    pcaColnames(tS)
    hclustColnames(tS)

    tS = timeSeries(matrix(rnorm(500), 100), units = sample(LETTERS[1:5]),
        format = "counts")
    head(tS)
    colnames(tS)
    sortColnames(tS)
    sampleColnames(tS)
    orderColnames(tS)
    statsColnames(tS, FUN = colMeans)
    pcaColnames(tS)
    hclustColnames(tS)
}


################################################################################

