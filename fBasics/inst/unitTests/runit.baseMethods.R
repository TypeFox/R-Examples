
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
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# FUNCTION:                 BASIC EXTENSIONS:                                         
#  align.default             align default method                                         
#  atoms.default             atoms default method                  
#  attach.default            attach default method           
#  colnames<-.default        colnames<- default method                  
#  cor.default               cor default method                      
#  cov.default               var default method                                          
#  log.default               log default method                  
#  outlier.default           outlier default method             
#  rownames<-.default        rownames<- default method
#  rank.default              rank default method                                                                     
#  sample.default            sample default method                                              
#  sort.default              sort default method                  
#  stdev.default             stdev default method                 
#  termPlot.default          termPlot default method                                               
#  var.default               var default method                  
#  volatility.default        volatility default method                                                                                                                                                                                                         
################################################################################


test.as.align.default <- 
    function() 
{
    NA
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.colnames.default = 
    function()
{
    #  "colnames<-.default"          colnames<- default method

    # Row and Column Names:
    m = matrix(1:8, ncol = 2)
    m
    
    # Set Names:
    colnames(m) <- c("A", "B")
    m
    
    # Get Names:
    colnames(m)
    
    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.rownames.default = 
    function()
{
    #  "rownames<-.default"          rownames<- default methodd

    # Row and Column Names:
    m = matrix(1:8, ncol = 2)
    m
    
    # Set Names:
    rownames(m) = as.character(1:4)
    m
    
    # Get Names:
    rownames(m)
    
    # Return Value:
    return()
}


################################################################################

