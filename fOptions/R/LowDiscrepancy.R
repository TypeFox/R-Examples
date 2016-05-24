
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
# FUNCTION:             DESCRIPTION:
#  runif.pseudo           Uniform Pseudo Random number sequence
#  rnorm.pseudo           Normal Pseudo Random number sequence
#  runif.halton           Uniform Halton low discrepancy sequence
#  rnorm.halton           Normal Halton low discrepancy sequence
#  runif.sobol            Uniform Sobol low discrepancy sequence
#  rnorm.sobol            Normal Sobol low discrepancy sequence
################################################################################


runif.pseudo <- 
  function(n, dimension, init = NULL) 
{
    # Description:
    #   Uniform Pseudo Random number sequence

    matrix(runif(n*dimension), ncol = dimension)
}


# ------------------------------------------------------------------------------


rnorm.pseudo <- 
  function(n, dimension, init = TRUE) 
{

    # Description:
    #   Normal Pseudo Random number sequence

    matrix(rnorm(n*dimension), ncol = dimension)
}


# -----------------------------------------------------------------------------


runif.halton <- 
  function (n, dimension, init = TRUE)
{   
    # A function implemented by Diethelm Wuertz

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
        .setfOptionsEnv(.runif.halton.seed = list(base = rep(0, dimension), offset = 0))
    }
    optEnv <- .getfOptionsEnv(".runif.halton.seed")

    # SUBROUTINE HALTON(QN, N, DIMEN, BASE, OFFSET, INIT, TRANSFORM)
    result <- .Fortran("halton",
                       qn = numeric(n*dimension),
                       as.integer(n),
                       as.integer(dimension),
                       base = as.integer(optEnv$base),
                       offset=as.integer(optEnv$offset),
                       as.integer(init),
                       0L,
                       PACKAGE = "fOptions")

    # For the next numbers save:
    .setfOptionsEnv(.runif.halton.seed = result[c("base", "offset")])

    matrix(result[["qn"]], ncol = dimension)
}


# ------------------------------------------------------------------------------


rnorm.halton <- function (n, dimension, init = TRUE)
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
        .setfOptionsEnv(.rnorm.halton.seed = list(base = rep(0, dimension), offset = 0))
    }
    optEnv <- .getfOptionsEnv(".rnorm.halton.seed")

    # SUBROUTINE HALTON(QN, N, DIMEN, BASE, OFFSET, INIT, TRANSFORM)
    result <- .Fortran("halton",
                       qn = numeric(n * dimension),
                       as.integer(n),
                       as.integer(dimension),
                       base = as.integer(optEnv$base),
                       offset=as.integer(optEnv$offset),
                       as.integer(init),
                       1L,
                       PACKAGE = "fOptions")

    # For the next numbers save:
    .setfOptionsEnv(.rnorm.halton.seed = result[c("base", "offset")])

    matrix(result[["qn"]], ncol = dimension)
}


# -----------------------------------------------------------------------------


runif.sobol <- function (n, dimension, init = TRUE, scrambling = 0, seed = 4711)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Uniform Sobol Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #           N : LD numbers to create
    #  SCRAMBLING : One of the numbers 0,1,2,3
    #

    # FUNCTION:
    stopifnot(0 <= (scrambling <- as.integer(scrambling)), scrambling <= 3)

    # Restart Settings:
    if (init) {
        .setfOptionsEnv(.runif.sobol.seed = list(quasi = rep(0, dimension), ll = 0,
                        count = 0, sv = rep(0, dimension*30), seed = seed))
    }
    optEnv <- .getfOptionsEnv(".runif.sobol.seed")

    # SSOBOL(QN,N,DIMEN,QUASI,LL,COUNT,SV,scrambling,SEED,INIT,TRANSFORM)
    result <- .Fortran("sobol",
                       qn = numeric(n * dimension),
                       as.integer(n),
                       as.integer(dimension),
                       quasi = as.double (optEnv$quasi),
                       ll    = as.integer(optEnv$ll),
                       count = as.integer(optEnv$count),
                       sv    = as.integer(optEnv$sv),
                       scrambling,
                       seed  = as.integer(optEnv$seed),
                       as.integer(init),
                       0L,
                       PACKAGE = "fOptions")

    # For the next numbers save:
    .setfOptionsEnv(.runif.sobol.seed = result[c("quasi","ll","count","sv","seed")])

    matrix(result[["qn"]], ncol = dimension)
}


# ------------------------------------------------------------------------------


rnorm.sobol <- function (n, dimension, init = TRUE, scrambling = 0, seed = 4711)
{   # A function implemented by Diethelm Wuertz

    # Description:
    #   Normal Sobol Low Discrepancy Sequence

    # Details:
    #   DIMENSION : dimension <= 200
    #           N : LD numbers to create
    #  SCRAMBLING : One of the numbers 0,1,2,3

    # FUNCTION:

    stopifnot(0 <= (scrambling <- as.integer(scrambling)), scrambling <= 3)
    # Restart Settings:
    if (init) {
        .setfOptionsEnv(.rnorm.sobol.seed = list( quasi = rep(0, dimension), ll = 0,
                        count = 0, sv = rep(0, dimension*30), seed = seed))
    }
    optEnv <- .getfOptionsEnv(".rnorm.sobol.seed")

    # SSOBOL(QN,N,DIMEN,QUASI,LL,COUNT,SV,scrambling,SEED,INIT,TRANSFORM)
    result <- .Fortran("sobol",
                       qn = numeric(n * dimension),
                       as.integer(n),
                       as.integer(dimension),
                       quasi = as.double (optEnv$quasi),
                       ll    = as.integer(optEnv$ll),
                       count = as.integer(optEnv$count),
                       sv    = as.integer(optEnv$sv),
                       scrambling,
                       seed  = as.integer(optEnv$seed),
                       as.integer(init),
                       1L,
                       PACKAGE = "fOptions")

    # For the next numbers save:
    .setfOptionsEnv(.rnorm.sobol.seed = result[c("quasi","ll","count","sv","seed")])

    matrix(result[["qn"]], ncol = dimension)
}


################################################################################

