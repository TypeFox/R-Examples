
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
#   1999 - 2008, Diethelm Wuertz, Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               SOLVER:
#  .garchRnlminb           R coded solver nlmin
#  .garchRlbgfsb           R coded solver optim using method lbgfsb
#  .garchRnm               R coded solver nm as hybrid addon
#  .garchFsqp              Fortran coded solver sqp
#  .garchRdonlp2           R coded solver donlp2
#  .garchFmnfb             Fortran coded solver mnfb
################################################################################


## test.garchSolver.dem2gbp <-
## function()
## {
##     # Note, Default has changed: "cda" -> "fda"

##     # Loda Data
##     data(dem2gbp)

##     garchFit(~ garch(1,1), dem2gbp, hessian = "fda", algorithm = "nlminb")
##     # Time difference of 2.969 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006190    0.008464   -0.731 0.464545
##     # omega   0.010761    0.003081    3.492 0.000479 ***
##     # alpha1  0.153134    0.027720    5.524 3.31e-08 ***
##     # beta1   0.805974    0.035761   22.537  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 *** 0.001 ** 0.01 * 0.05 . 0.1   1
##     #
##     # Log Likelihood:
##     #  1106.608    normalized:  0.5605916

##     garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "nlminb")

##     garchFit(~ garch(1,1), dem2gbp, algorithm = "lbfgsb")
##     # Time difference of 4.75 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006187    0.008462   -0.731 0.464720
##     # omega   0.010785    0.002860    3.771 0.000163 ***
##     # alpha1  0.153306    0.026556    5.773 7.79e-09 ***
##     # beta1   0.805706    0.033612   23.971  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  1106.608    normalized:  0.5605916


##     garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "nlminb+nm")
##     # Time difference of 6.094 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006190    0.008462   -0.732 0.464447
##     # omega   0.010761    0.002853    3.772 0.000162 ***
##     # alpha1  0.153134    0.026523    5.774 7.76e-09 ***
##     # beta1   0.805974    0.033553   24.021  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  1106.608    normalized:  0.5605916
##     #


##     garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "lbfgsb+nm")


##     garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "donlp2")
##     # Time difference of 7.094 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006190    0.008462   -0.732 0.464447
##     # omega   0.010761    0.002853    3.772 0.000162 ***
##     # alpha1  0.153134    0.026523    5.774 7.76e-09 ***
##     # beta1   0.805973    0.033553   24.021  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  1106.608    normalized:  0.5605916
##     #


## garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "sqp")
##     # Time difference of 1.906 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006190    0.008462   -0.732 0.464437
##     # omega   0.010761    0.002853    3.772 0.000162 ***
##     # alpha1  0.153134    0.026523    5.774 7.76e-09 ***
##     # beta1   0.805974    0.033553   24.021  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  0.5605916    normalized:  0.0002839877
##     #

##     garchFit(~ garch(1,1), dem2gbp, hessian = "cda", algorithm = "mnfb")
##     # Time difference of 1.344 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu     -0.006190    0.008461   -0.732 0.464384
##     # omega   0.010761    0.002853    3.772 0.000162 ***
##     # alpha1  0.153134    0.026523    5.774 7.75e-09 ***
##     # beta1   0.805974    0.033552   24.021  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  0.5605916    normalized:  0.0002839877

##     # Return Value:
##     return()
## }

# ------------------------------------------------------------------------------

test.garchSolver2.dem2gbp <-
function()
{
    # Note, Default has changed: "cda" -> "fda"

    # Loda Data
    data(dem2gbp)

##    garchFit(~ garch(1,1), dem2gbp, hessian = "fcd", algorithm = "nlminb")

    garchFit(~ garch(1,1), dem2gbp, hessian = "rcd", algorithm = "nlminb")

    garchFit(~ garch(1,1), dem2gbp, algorithm = "lbfgsb")

##    garchFit(~ garch(1,1), dem2gbp, hessian = "fcd", algorithm = "nlminb+nm")

    garchFit(~ garch(1,1), dem2gbp, hessian = "rcd", algorithm = "lbfgsb+nm")

    # garchFit(~ garch(1,1), dem2gbp, hessian = "fcd", algorithm = "donlp2")

#    garchFit(~ garch(1,1), dem2gbp, hessian = "rcd", algorithm = "sqp")

  ##  garchFit(~ garch(1,1), dem2gbp, hessian = "fcd", algorithm = "mnfb")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


## test.garchSolver.sp500dge <-
## function()
## {
##     # Loda Data:
##     data(sp500dge)
##     sp500dge = 100 * sp500dge


##     garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "cda",
##         algorithm = "lbfgsb")


   #  garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "cda",
   #      algorithm = "sqp")
##     # Time difference of 48.954 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu      0.020646    0.006350    3.251  0.00115 **
##     # ma1     0.144745    0.008358   17.319  < 2e-16 ***
##     # omega   0.009988    0.001085    9.201  < 2e-16 ***
##     # alpha1  0.083803    0.004471   18.742  < 2e-16 ***
##     # gamma1  0.373092    0.027997   13.326  < 2e-16 ***
##     # beta1   0.919401    0.004093  224.614  < 2e-16 ***
##     # delta   1.435124    0.067203   21.355  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     # 1.264345    normalized:  7.413338e-05
##     #

    # garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "fda",
    #     algorithm = "sqp")


    # garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "cda",
    #     algorithm = "mnfb")
##     # Time difference of 35.156 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu      0.020646    0.006351    3.251  0.00115 **
##     # ma1     0.144745    0.008358   17.319  < 2e-16 ***
##     # omega   0.009988    0.001086    9.201  < 2e-16 ***
##     # alpha1  0.083803    0.004472   18.741  < 2e-16 ***
##     # gamma1  0.373092    0.027997   13.326  < 2e-16 ***
##     # beta1   0.919401    0.004093  224.606  < 2e-16 ***
##     # delta   1.435124    0.067204   21.355  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  1.264345    normalized:  7.413338e-05
##     #


    # garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "fda",
    #     algorithm = "mnfb")
##     # Time difference of 31.516 secs
##     # Error Analysis:
##     #         Estimate  Std. Error  t value Pr(>|t|)
##     # mu      0.020646    0.006342    3.255  0.00113 **
##     # ma1     0.144745    0.008363   17.307  < 2e-16 ***
##     # omega   0.009988    0.001113    8.975  < 2e-16 ***
##     # alpha1  0.083803    0.004534   18.485  < 2e-16 ***
##     # gamma1  0.373092    0.027798   13.422  < 2e-16 ***
##     # beta1   0.919401    0.004102  224.141  < 2e-16 ***
##     # delta   1.435124    0.066283   21.652  < 2e-16 ***
##     # ---
##     # Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
##     #
##     # Log Likelihood:
##     #  1.264345    normalized:  7.413338e-05
##     #


##     # Return Value:
##     return()
## }

# ------------------------------------------------------------------------------

test.garchSolver.sp500dge <-
function()
{
    # Loda Data:
    data(sp500dge)
    sp500dge = 100 * sp500dge


#    garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "fcd",
#             algorithm = "lbfgsb")


###     garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "fcd",
###              algorithm = "sqp")

###     garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "rcd",
###              algorithm = "sqp")


#    garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "fcd",
#             algorithm = "mnfb")


#    garchFit(~ arma(0,1) + aparch(1,1), sp500dge, hessian = "rcd",
#             algorithm = "mnfb")


    garchFit(~ arma(0,1) + aparch(1,1), sp500dge)

    # Return Value:
    return()
}

################################################################################
