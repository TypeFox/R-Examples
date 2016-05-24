
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR Description. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA 02111-1307 USA


################################################################################
# FUNCTION:                    DESCRIPTION:
#  setWindowsFun<-              Sets name of rolling windows function
#  setWindowsParams<-           Sets additional parameters to windows function
#  setWindowsHorizon<-          Sets horizon of the rolling window
# FUNCTION:                    DESCRIPTION:         
#  setStrategyFun<-             Sets name of portfolio strategy function    
#  setStrategyParams<-          Sets additional parameters to strategy function
# FUNCTION:                    DESCRIPTION:
#  setSmootherFun<-             Sets name of weights smoothing function
#  setSmootherParams<-          Sets additional parameters to smoother function
#  setSmootherLambda<-          Sets lambda for EMA smoothing
#  setSmootherDoubleSmoothing<- Sets double ema setting, logical
#  setSmootherInitialWeights<-  Sets initial weights of the portfolio
#  setSmootherSkip<-            Sets number of months to skip starting
################################################################################


"setWindowsFun<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets name of rolling windows function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@windows$windows <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setWindowsParams<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets additional parameters to windows function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@windows$params <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setWindowsHorizon<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets horizon of the rolling window
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@windows$params$horizon <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setStrategyFun<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets portfolio strategy function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@strategy$strategy <- value
    
    # Return Value:
    backtest
}

# ------------------------------------------------------------------------------
          

"setStrategyParams<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets additional parameters to strategy function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@strategy$params <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherFun<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets name of weights smoothing function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$smoother <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherParams<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets additional parameters to smoother function
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$params <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherLambda<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets lambda parameter for EMA smoothing
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$params$lambda <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherDoubleSmoothing<-" <- 
function(backtest, value)   
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets double EMA setting, TRUE or FALSE, a logical
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$params$doubleSmoothing <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherInitialWeights<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets initial weights of the portfolio
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$params$initialWeights <- value
    
    # Return Value:
    backtest
}


# ------------------------------------------------------------------------------


"setSmootherSkip<-" <- 
function(backtest, value)
{
    # A function implemented by William Chen
    
    # Description:
    #   Sets number of months to skip starting values
    
    # Arguments:
    
    # FUNCTION:
    
    # Set Value:
    backtest@smoother$params$skip <- value
    
    # Return Value:
    backtest
}
   

################################################################################
  
 