
# Copyrights (C) 2014 Thiago do Rego Sousa <thiagoestatistico@gmail.com>

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
#  FUNCTION:               DESCRIPTION:
#
#  .getOrder        Return a matrix with parameter order to be used 
#                          inside function GSgarch.FitAIC 
################################################################################


.getOrder <- 
    function(order.max = c(1,1,1,1))
{
      
    # Description:
    #   Iterates over the parameters to create a vector with parameter
    #   orders (like (1,1,0,1)) to use inside function GSgarch.FitAIC
    
    # Arguments:
    #   nMAX, mMAX, pMAX, qMAX - maximum order to be estimated
    
    # Return:
    #   arma.garch.order - A matrix with several model orders in the 
    #   format [m,n,p,q] 
    
    # FUNCTION:            
      
    # error treatment on input parameters
    m = order.max[1]; n = order.max[2]; p = order.max[3]; q = order.max[4] 
    if(m%%1 != 0 || n%%1 != 0 || p%%1 != 0 || q%%1 != 0 || 
        any (c(m,n,p,q) < 0) || (p == 0 && q != 0) ||
        any (c(m,n,p,q) > 10) ) 
        stop ("Invalid ARMA-GARCH order. We allow pure GARCH or APARCH. AR/MA/ARMA-GARCH/APARCH models.
              The order of the parameters could be set up to 10.")
    arma.garch.order <- c()
    for(i1 in 0:m)
    {
        for(i2 in 0:n)
        {
            for(i3 in 1:p)
            {
                for(i4 in 0:q)
                { 
                    ord <- c(i1,i2,i3,i4)
                    arma.garch.order <- rbind(arma.garch.order,c(i1,i2,i3,i4))
                }
            }
        }
    }
    return(arma.garch.order)
}


################################################################################

