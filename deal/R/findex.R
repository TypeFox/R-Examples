## findex.R
## Author          : Claus Dethlefsen
## Created On      : Thu Nov 29 10:15:11 2001
## Last Modified By: Claus Dethlefsen
## Last Modified On: Tue Jul 22 16:53:14 2003
## Update Count    : 67
## Status          : Unknown, Use with caution!
###############################################################################
##
##    Copyright (C) 2002  Susanne Gammelgaard Bøttcher, Claus Dethlefsen
##
##    This program is free software; you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation; either version 2 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program; if not, write to the Free Software
##    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
######################################################################

findex <- function(i, dim, config=TRUE) {
    ## find index for use with an array of dimension 'dim'
    ##
    ## if config==T :: (i is a configuration matrix)
    ## i is then interpreted as a
    ## matrix with one row per wanted entry. The columns are the
    ## configurations of each of the discrete variables (in the proper
    ## order).
    ## Returned is a vector of length the number of rows of i. The
    ## entries correspond to each row and is the corresponding number if
    ## the array were 'folded' out.
    ##
    ## if config==F ::
    ## i is a vector of indices in the unfolded array. We want the
    ## corresponding configurations of the discrete variables
    ## output is a matrix with one row per configuration
    ## 
    ## Thus, findex(config=T) and findex(config=F) are each others
    ## inverse functions
    
    mymod <- function(a,n) ifelse(a%%n==0,a%%n+n,a%%n)
    
    roundup <- function(a) floor(a+0.999)
    
    
    N <- prod(dim)
    D <- length(dim)
    
    if (config) res <- array(1:N,dim=dim)[i]
    
    else {
        ## Like V&R page 42
        
        res <- matrix(NA,length(i),D)
        for (k in 1:length(i)) {
            j <- i[k]
            res[k,1] <- mymod(j,dim[1])
            if (D>1) { 
                for (s in 2:D) 
                    res[k,s] <-  roundup(mymod(j,prod(dim[1:s]))/
                                         prod(dim[1:(s-1)]))
            }
        }
    }
    
    res
}

