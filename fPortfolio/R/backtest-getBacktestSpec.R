
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
# FUNCTION:                    DESCRIPTION:
#  getWindows                   Extracts windows slot
#   getWindowsFun               Extracts name of windows function
#   getWindowsParams            Extracts a list of windows specific parameters
#   getWindowsHorizon           Extracts windows horizon
# FUNCTION:                    DESCRIPTION:
#  getStrategy                  Extracts strategy slot
#   getStrategyFun              Extracts the name of portfolio strategy function
#   getStrategyParams           Extracts a list of strategy specific parameters
# FUNCTION:                    DESCRIPTION:
#  getSmoother                  Extracts the smoother slot
#   getSmootherFun              Extracts the name of the moother function
#   getSmootherParams           Extracts a list of smoothing specific parameters
#   getSmootherLambda           Extracts the smoothing parameter Lambda
#   getSmootherDoubleSmoothing  Extracts setting for double smoothing
#   getSmootherInitialWeights   Extracts the initial weights in the smoothing
#   getSmootherSkip             Extracts the number of skipped months
# FUNCTION:                    DESCRIPTION:
#  getMessages                  Extracts the message slot
################################################################################


getWindows.fPFOLIOBACKTEST <- 
    function(object)
{   
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts Windows slot from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
        
    # Description:
    #   gets the "model" slot from an object of class 4
   
    # Arguments:
    #   object - an object of class S4
   
    # FUNCTION:
   
    # Return Value:
    getSlot(object, "windows")
}
            

# ------------------------------------------------------------------------------
    

getWindowsFun.fPFOLIOBACKTEST <- 
function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts name of windows function from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getWindows(object)$windows
}


# ------------------------------------------------------------------------------


getWindowsParams.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts a list of windows specific parameters from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getWindows(object)$params
}

   
# ------------------------------------------------------------------------------


getWindowsHorizon.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts windows horizon from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getWindowsParams(object)$horizon
}
    
    
# ------------------------------------------------------------------------------


getSmoother.fPFOLIOBACKTEST <- 
    function(object)
{   
    # A function implemented by William Chen and Diethelm Wuertz
        
    # Description:
    #   Extracts the smoother slot from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
   
    # Return Value:
    getSlot(object, "smoother")
}


# ------------------------------------------------------------------------------


getSmootherFun.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the name of the moother function from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmoother(object)$smoother
}

    
# ------------------------------------------------------------------------------


getSmootherParams.fPFOLIOBACKTEST <-  
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts a list of strategy specific parameters
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmoother(object)$params
}

    
# ------------------------------------------------------------------------------


getSmootherLambda.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the smoothing parameter Lambda from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmootherParams(object)$lambda
}

    
# ------------------------------------------------------------------------------


getSmootherDoubleSmoothing.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts setting for double smoothing from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmootherParams(object)$doubleSmoothing
}

    
# ------------------------------------------------------------------------------


getSmootherInitialWeights.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the initial weights in the smoothing from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmootherParams(object)$initialWeights
}

    
# ------------------------------------------------------------------------------


getSmootherSkip.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the number of skipped months from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getSmootherParams(object)$skip
}
 
       
# ------------------------------------------------------------------------------


getStrategy.fPFOLIOBACKTEST <- 
    function(object)
{   
    # A function implemented by William Chen and Diethelm Wuertz
        
    # Description:
    #   Extracts strategy slot from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
   
    # Return Value:
    getSlot(object, "strategy")
}

    
# ------------------------------------------------------------------------------


getStrategyFun.fPFOLIOBACKTEST <-  
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the name of portfolio strategy function from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getStrategy(object)$strategy
}

    
# ------------------------------------------------------------------------------


getStrategyParams.fPFOLIOBACKTEST <- 
    function(object) 
{
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts a list of strategy specific parameters from an object 
    #   of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
    
    # Return Value:
    getStrategy(object)$params
}
   

# ------------------------------------------------------------------------------
                

getMessages.fPFOLIOBACKTEST <- 
    function(object)
{   
    # A function implemented by William Chen and Diethelm Wuertz
    
    # Description:
    #   Extracts the message slot from an object of class fPFOLIOBACKTEST
    
    # Arguments:
    #   object - an object of class fPFOLIOBACKTEST
    
    # FUNCTION:
   
    # Return Value:
    getSlot(object, "messages")
}


################################################################################

