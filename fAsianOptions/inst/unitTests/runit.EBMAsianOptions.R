
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
# MOMENT MATCHING:                   DESCRIPTION:
#  MomentMatchedAsianOption           Valuate moment matched option prices
#  .LevyTurnbullWakemanAsianOption     Log-Normal Approximation
#  .MilevskyPosnerAsianOption          Reciprocal-Gamma Approximation
#  .PosnerMilevskyAsianOption          Johnson Type I Approximation
#  MomentMatchedAsianDensity          Valuate moment matched option densities
#  .LevyTurnbullWakemanAsianDensity    Log-Normal Approximation
#  .MilevskyPosnerAsianDensity         Reciprocal-Gamma Approximation
#  .PosnerMilevskyAsianDensity         Johnson Type I Approximation 
# GRAM CHARLIER SERIES EXPANSION:    DESCRIPTION:
#  GramCharlierAsianOption            Calculate Gram-Charlier option prices
#  .GramCharlierAsianDensity          NA
# STATE SPACE MOMENTS:               DESCRIPTION:
#  AsianOptionMoments                 Methods to calculate Asian Moments
#  .DufresneAsianOptionMoments         Moments from Dufresne's Formula
#  .AbrahamsonAsianOptionMoments       Moments from Abrahamson's Formula
#  .TurnbullWakemanAsianOptionMoments  First 2 Moments from Turnbull-Wakeman
#  .TolmatzAsianOptionMoments          Asymptotic Behavior after Tolmatz
# STATE SPACE DENSITIES:              DESCRIPTION:
#  StateSpaceAsianDensity              NA
#  .Schroeder1AsianDensity             NA
#  .Schroeder2AsianDensity             NA
#  .Yor1AsianDensity                   NA
#  .Yor2AsianDensity                   NA
#  .TolmatzAsianDensity                NA
#  .TolmatzAsianProbability            NA
# PARTIAL DIFFERENTIAL EQUATIONS:     DESCRIPTION:
#  PDEAsianOption                      PDE Asian Option Pricing
#   .ZhangAsianOption                   Asian option price by Zhang's 1D PDE
#    ZhangApproximateAsianOption
#   .VecerAsianOption                   Asian option price by Vecer's 1D PDE 
# LAPLACE INVERSION:                  DESCRIPTION:
#   GemanYorAsianOption                Asian option price by Laplace Inversion
#   gGemanYor                          Function to be Laplace inverted
# SPECTRAL EXPANSION:                 DESCRIPTION:
#   LinetzkyAsianOption                Asian option price by Spectral Expansion
#   gLinetzky                          Function to be integrated
# BOUNDS ON OPTION PRICES:            DESCRIPTION:
#   BoundsOnAsianOption                 Lower and upper bonds on Asian calls
#    CurranThompsonAsianOption          From Thompson's continuous limit
#    RogerShiThompsonAsianOption        From Thompson's single integral formula   
#    ThompsonAsianOption                Thompson's upper bound 
# SYMMETRY RELATIONS:                 DESCRIPTION:
#   CallPutParityAsianOption           Call-Put parity Relation
#   WithDividendsAsianOption           Adds dividends to Asian Option Formula
# TABULATED RESULTS:                  DESCRIPTION:
#   FuMadanWangTable                   Table from Fu, Madan and Wang's paper
#   FusaiTaglianiTable                 Table from Fusai und tagliani's paper
#   GemanTable                         Table from Geman's paper
#   LinetzkyTable                      Table from Linetzky's paper
#   ZhangTable                         Table from Zhang's paper
#   ZhangLongTable                     Long Table from Zhang's paper
#   ZhangShortTable                    Short Table from Zhang's paper
################################################################################


test.MomentMatchedAsianOption = 
function()
{
    # MomentMatchedAsianOption

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.MomentMatchedAsianDensity = 
function()
{
    # MomentMatchedAsianDensity  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.GramCharlierAsianOption = 
function()
{
    # GramCharlierAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.AsianOptionMoments = 
function()
{
    # AsianOptionMoments  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.StateSpaceAsianDensity = 
function()
{
    # StateSpaceAsianDensity  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.PDEAsianOption = 
function()
{
    # PDEAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.GemanYorAsianOption = 
function()
{
    # GemanYorAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.LinetzkyAsianOption = 
function()
{
    # LinetzkyAsianOption  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.BoundsOnAsianOption = 
function()
{    
    # BoundsOnAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.CurranThompsonAsianOption = 
function()
{
    # CurranThompsonAsianOption  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.RogerShiThompsonAsianOption = 
function()
{
    # RogerShiThompsonAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ThompsonAsianOption = 
function()
{
    # ThompsonAsianOption
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.CallPutParityAsianOption = 
function()
{
    # CallPutParityAsianOption  

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.WithDividendsAsianOption = 
function()
{
    # WithDividendsAsianOption 

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.Table = 
function()
{ 
    # Table 

    # Return Value:
    return()    
}

    
# ------------------------------------------------------------------------------


if (FALSE) {
    require(RUnit)
    testResult <- runTestFile("C:/Rmetrics/SVN/trunk/fOptions/tests/runit3E.R",
        rngKind = "Marsaglia-Multicarry", rngNormalKind = "Inversion")
    printTextProtocol(testResult)
}


################################################################################

