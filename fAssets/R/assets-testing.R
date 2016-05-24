
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                   DESCRIPTION:
#  assetsTest                  Tests for multivariate Normal Assets
# FUNCTION:                   DESCRIPTION:
#  mvshapiroTest               Multivariate Shapiro Test
#  mvenergyTest                Multivariate E-Statistic (Energy) Test
################################################################################


assetsTest <-
    function(x, method = c("shapiro", "energy"), Replicates = 99)
{
    # Description:
    #   Tests for multivariate Normal Assets

    # Example:
    #   mvnormTest(x = assetsSim(100))
    #   mvnormTest(x = assetsSim(100), method = "e", Replicates = 99)

    # FUNCTION:

    # Test:
    method <- match.arg(method)
    if (method == "shapiro") {
        test <- mvshapiroTest(x)
    }
    if (method == "energy") {
        test <- mvenergyTest(x, Replicates = Replicates)
    }

    # Return Value:
    test
}


# ------------------------------------------------------------------------------


mvshapiroTest <-
    function(x)
{
    # Description:
    #   Computes Shapiro's normality test for multivariate variables
    
    # Requires: 
    #   Package: mvnormtest
    #   Version: 0.1-6
    #   Date: 2005-04-02
    #   Title: Normality test for multivariate variables
    #   Author: Slawomir Jarek
    #   Maintainer: Slawomir Jarek <slawomir.jarek@gallus.edu.pl>
    #   Description: Generalization of shapiro-wilk test for
    #       multivariate variables.
    #   License: GPL
              
    # Example:
    #   mvshapiroTest(x = assetsSim(100))

    # FUNCTION:

    # Transform:
    U <- t(as.matrix(x))
    
    # Test
    test <- mvnormtest::mshapiro.test(U)
    
    # Return Value:
    test
}


# ------------------------------------------------------------------------------


mvenergyTest <- 
    function(x, Replicates = 99)
{
    # Description:
    #   Computes E-statistics test for multivariate variables

    # Requires:
    #   Package: energy
    #   Author: Maria L. Rizzo <mrizzo @ bgnet.bgsu.edu> and
    #       Gabor J. Szekely <gabors @ bgnet.bgsu.edu>
    #   License: GPL 2.0 or later

    # Example:
    #   mvenergyTest(x = assetsSim(100), 99)

    # FUNCTION:

    # Transform:
    x <- as.matrix(x)
    
    # Test:
    test <- energy::mvnorm.etest(x, R = Replicates)
    
    # Return Value:
    test
}


################################################################################

