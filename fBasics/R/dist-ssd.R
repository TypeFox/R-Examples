
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
# FUNCTIONS:       DESCRIPTION:
#  dssd             Returns spline smoothed density estimate
#  pssd             Returns spline smoothed probability estimate
#  qssd             Returns spline smoothed quantiles estimate
#  rssd             Returns spline smoothed random variates 
# LIBRARY:             DESCRPTION:
#  stabledist           Stable Duistribution
################################################################################


dssd <- 
   function(x, param, log = FALSE) 
{    
    # A function implemented by Diethelm Wuertz

    # Description:
    #  Evaluate density using smoothing spline ANOVA model 
    
    # FUNCTION:
    
    # Density:
    ans <- gss::dssden(object = param, x = x) 
    
    # Log:
    if(log) ans <- log(ans)
    
    # Return Value:
    ans
}    


# ------------------------------------------------------------------------------

 
pssd <- 
    function(q, param) 
{    
    # A function implemented by Diethelm Wuertz

    # Description:
    #      Evaluate probability using smoothing spline ANOVA model 
    
    # FUNCTION:
    
    # Return Value:
    gss::pssden(object = param, q = q) 
}


# ------------------------------------------------------------------------------


qssd <- 
    function(p, param) 
{    
    # A function implemented by Diethelm Wuertz

    # Description:
    #     Evaluate quantiles using smoothing spline ANOVA model     
    
    # FUNCTION:
    
    # Return Value:
    gss::qssden(object = param, p = p) 
}


# ------------------------------------------------------------------------------


rssd <-  
    function(n, param) 
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #     Generate random deviates using smoothing spline ANOVA model 
    
    # FUNCTION:
    
    # Return Value:
    gss::qssden(object = param, p = runif(n)) 
}


################################################################################


