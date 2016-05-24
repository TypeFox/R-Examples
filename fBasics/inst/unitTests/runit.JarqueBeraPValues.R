
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
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:           JARQUE-BERA DATA TABLES:
# .jbTable             Finite sample p values for the Jarque Bera test
# .jbPlot              Plots probability
# .pjb                 Returns probabilities for the JB Test given quantiles
# .qjb                 Returns quantiles for the ADF Test given probabilities
# DATA:               Description:
# .jbLM                Jarque-Bera Lagrange Multiplier Test Data
# .jbALM               Jarque Bera Augmented Lagrange Multiplier Test Data
################################################################################


test.jbTable =
function()
{
    if (FALSE) {

        require(akima)

        # Jarque-Bera Table:
        #   .jbTable(type = c("LM", "ALM"), size = c("mini", "small", "all"))
        table = .jbTable()
        table

        # Perspective Plot:
        #   .jbPlot(type = c("LM", "ALM"))
        .jbPlot()

    }

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.jbData =
function()
{
    # Jarque-Bera LM Data:
    class(fBasics:::.jbLM)
    class(fBasics:::.jbLM())
    head(fBasics:::.jbLM())

    # Jarque-Bera ALM Data:
    class(fBasics:::.jbALM)
    class(fBasics:::.jbALM())
    head(fBasics:::.jbALM())

    # Return Value:
    return()
}


################################################################################

