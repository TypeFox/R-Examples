# LatticeKrig  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2012
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2

LKDistGridComponents<- function( x1, gridList, delta,
                         max.points = NULL, 
                      mean.neighbor = NULL,
                      distance.type = "Euclidean"
                     ){
         LKGridCheck( distance.type, x1, gridList) 
         n1<- nrow(x1)
         info<- summary( gridList ) 
# figure out how large an array to allocate for distance matrix 
         Nmax<- LKGridFindNmax(n1, max.points, mean.neighbor,delta,gridList)        
# NOTE all the dx's should be equal -- this is checked above 
# FORTRAN expects x1 to be centered and scaled to the grid. 
# i.e. the grid coordinates  are all just integer values 1:10 etc.
         xScaled <- scale( x1, center = info$min, scale = info$dx) + 1 
         deltaScaled <- delta/mean(info$dx)    	
    out <- .Fortran("LKdistgridcomp",
                                  x1 = as.double(xScaled),
                                  n1 = as.integer(n1), 
                               nGrid = as.integer( info$n),
                                nDim = as.integer(info$dim), 
                               delta = as.double(deltaScaled),
                                irow = as.integer(rep(0, Nmax)),
                                jcol = as.integer(rep(0, Nmax)),       
                                  ra = as.double(rep(-1, Nmax*info$dim)),
                                Nmax = as.integer(Nmax),
                               iflag = as.integer(1),
                               index = as.integer( rep(-1,Nmax*info$dim)),
                             PACKAGE = "LatticeKrig")
                             
# out$Nmax are now the actual number of nonzero distances found. 
# unless there was an error (iflag!=0)                            
      N <- out$Nmax
      if (N==0) {
            stop("delta too small! -- no points found")
            }
      if( out$iflag <0){
      	stop( paste( "Nmax = " , N, " exceeded")  )
      }           
      DX<- mean(info$dx)
      raOut<- DX* matrix(out$ra,  nrow=Nmax, ncol=info$dim )[1:N,]                     
  	  out<- list(
  	             ind = cbind( out$irow[1:N], out$jcol[1:N] ) ,
                  ra = raOut,                  
                  da = c(n1,prod(out$nGrid))
                )            
               
      return(out)
    }
 
