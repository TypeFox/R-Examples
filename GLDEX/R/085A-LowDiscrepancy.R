
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
#   1999 - 2004, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:             DESCRIPTION:
#  runif.pseudo           Uniform Pseudo Random number sequence
#  rnorm.pseudo           Normal Pseudo Random number sequence
#  runif.halton           Uniform Halton low discrepancy sequence
#  rnorm.halton           Normal Halton low discrepancy sequence
#  runif.sobol            Uniform Sobol low discrepancy sequence
#  rnorm.sobol            Normal Sobol low discrepancy sequence
################################################################################

runif.pseudo =
function(n, dimension, init = NULL)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Uniform Pseudo Random number sequence

    # FUNCTION:

    # Deviates:
    result = matrix(runif(n*dimension), ncol = dimension)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rnorm.pseudo =
function(n, dimension, init = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Normal Pseudo Random number sequence

    # FUNCTION:

    # Deviates:
    result = matrix(rnorm(n*dimension), ncol = dimension)

    # Return Value:
    result
}


# -----------------------------------------------------------------------------


runif.halton =
function (n, dimension, init = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Uniform Halton Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #       N : LD numbers to create

    # FUNCTION:

    # Restart Settings:
    if (init) {
        ## YC : this code should not needed since we have now global Env
        # .warn = options()$warn
        # options(warn = -1)
        # rm(".runif.halton.seed")
        # options(warn = .warn)
        .setGLDEXEnv(.runif.halton.seed = list(base = rep(0, dimension), offset = 0))
    }

    # Generate:
    qn = rep(0, n*dimension)

    # SUBROUTINE HALTON(QN, N, DIMEN, BASE, OFFSET, INIT, TRANSFORM)
    result = .Fortran("halton",
        as.double(qn),
        as.integer(n),
        as.integer(dimension),
        as.integer(.getGLDEXEnv(".runif.halton.seed")$base),
        as.integer(.getGLDEXEnv(".runif.halton.seed")$offset),
        as.integer(init),
        as.integer(0),
        PACKAGE = "GLDEX")

    # For the next numbers save:
###     .warn = options()$warn
###     options(warn = -1)
###     rm(".runif.halton.seed")
###     options(warn = .warn)
    .setGLDEXEnv(.runif.halton.seed = list(base = result[[4]], offset = result[[5]]))

    # Deviates:
    result = matrix(result[[1]], ncol = dimension)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rnorm.halton =
function (n, dimension, init = TRUE)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Normal Halton Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #       N : LD numbers to create

    # FUNCTION:

    # Restart Settings:
    if (init) {
###         .warn = options()$warn
###         options(warn = -1)
###         rm(".rnorm.halton.seed")
###         options(warn = .warn)
        .setGLDEXEnv(.rnorm.halton.seed = list(base = rep(0, dimension), offset = 0))
    }

    # Generate:
    qn = rep(0, n*dimension)

    # SUBROUTINE HALTON(QN, N, DIMEN, BASE, OFFSET, INIT, TRANSFORM)
    result = .Fortran("halton",
        as.double(qn),
        as.integer(n),
        as.integer(dimension),
        as.integer(.getGLDEXEnv(".rnorm.halton.seed")$base),
        as.integer(.getGLDEXEnv(".rnorm.halton.seed")$offset),
        as.integer(init),
        as.integer(1),
        PACKAGE = "GLDEX")

    # For the next numbers save:
###     .warn = options()$warn
###     options(warn = -1)
###     rm(".rnorm.halton.seed")
###     options(warn = .warn)
    .setGLDEXEnv(.rnorm.halton.seed = list(base = result[[4]], offset = result[[5]]))

    # Deviates:
    result = matrix(result[[1]], ncol = dimension)

    # Return Value:
    result
}


# -----------------------------------------------------------------------------


runif.sobol =
function (n, dimension, init = TRUE, scrambling = 0, seed = 4711)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Uniform Sobol Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #           N : LD numbers to create
    #  SCRAMBLING : One of the numbers 0,1,2,3
    #

    # FUNCTION:

    # Restart Settings:
    if (init) {
###         .warn = options()$warn
###         options(warn = -1)
###         rm(".runif.sobol.seed")
###         options(warn = .warn)
        .setGLDEXEnv(.runif.sobol.seed = list(quasi = rep(0, dimension), ll = 0,
            count = 0, sv = rep(0, dimension*30), seed = seed))
    }

    # Generate:
    qn = rep(0.0, n*dimension)

    # SSOBOL(QN,N,DIMEN,QUASI,LL,COUNT,SV,IFLAG,SEED,INIT,TRANSFORM)
    result = .Fortran("sobol",
        as.double(qn),
        as.integer(n),
        as.integer(dimension),
        as.double (.getGLDEXEnv(".runif.sobol.seed")$quasi),
        as.integer(.getGLDEXEnv(".runif.sobol.seed")$ll),
        as.integer(.getGLDEXEnv(".runif.sobol.seed")$count),
        as.integer(.getGLDEXEnv(".runif.sobol.seed")$sv),
        as.integer(scrambling),
        as.integer(.getGLDEXEnv(".runif.sobol.seed")$seed),
        as.integer(init),
        as.integer(0),
        PACKAGE = "GLDEX")

    # For the next numbers save:
###     .warn = options()$warn
###     options(warn = -1)
###     rm(".runif.sobol.seed")
###     options(warn = .warn)
    .setGLDEXEnv(.runif.sobol.seed = list(quasi = result[[4]], ll = result[[5]],
        count = result[[6]], sv = result[[7]], seed = result[[9]]))

    # Deviates:
    result = matrix(result[[1]], ncol = dimension)

    # Return Value:
    result
}


# ------------------------------------------------------------------------------


rnorm.sobol =
function (n, dimension, init = TRUE, scrambling = 0, seed = 4711)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Normal Sobol Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #           N : LD numbers to create
    #  SCRAMBLING : One of the numbers 0,1,2,3

    # FUNCTION:

    # Restart Settings:
    if (init) {
###         .warn = options()$warn
###         options(warn = -1)
###         rm(".rnorm.sobol.seed")
###         options(warn = .warn)
        .setGLDEXEnv(.rnorm.sobol.seed = list( quasi = rep(0, dimension), ll = 0,
            count = 0, sv = rep(0, dimension*30), seed = seed))
    }

    # Generate:
    qn = rep(0.0, n*dimension)

    # SSOBOL(QN,N,DIMEN,QUASI,LL,COUNT,SV,IFLAG,SEED,INIT,TRANSFORM)
    result = .Fortran("sobol",
        as.double(qn),
        as.integer(n),
        as.integer(dimension),
        as.double (.getGLDEXEnv(".rnorm.sobol.seed")$quasi),
        as.integer(.getGLDEXEnv(".rnorm.sobol.seed")$ll),
        as.integer(.getGLDEXEnv(".rnorm.sobol.seed")$count),
        as.integer(.getGLDEXEnv(".rnorm.sobol.seed")$sv),
        as.integer(scrambling),
        as.integer(.getGLDEXEnv(".rnorm.sobol.seed")$seed),
        as.integer(init),
        as.integer(1),
        PACKAGE = "GLDEX")

    # For the next numbers save:
###     .warn = options()$warn
###     options(warn = -1)
###     rm(".rnorm.sobol.seed")
###     options(warn = .warn)
    .setGLDEXEnv(.rnorm.sobol.seed = list(quasi = result[[4]], ll = result[[5]],
        count = result[[6]], sv = result[[7]], seed = result[[9]]))

    # Deviates:
    result = matrix(result[[1]], ncol = dimension)

    # Return Value:
    result
}


################################################################################

