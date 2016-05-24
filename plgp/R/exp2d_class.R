#******************************************************************************* 
#
# Particle Learning of Gaussian Processes
# Copyright (C) 2010, University of Cambridge
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (bobby@statslab.cam.ac.uk)
#
#*******************************************************************************

exp2d.C <- function(X, threed=TRUE)
  {
   if(is.null(X)) return(NULL);
   if(is.null(ncol(X))) X <- matrix(X, ncol=length(X))

   if(ncol(X) != 2)
     stop(paste("X should be a matrix (or data frame) with 2 cols, you have",
                ncol(X)))

   ## allocate space
   z <- rep(NA, nrow(X))
   
   ## for each row in X
   for(i in 1:nrow(X)) {
     
     ## extract ith coordinate
     x1 <- X[i,1]; x2 <- X[i,2]

     
     ## Hessian calculation
     E <- exp(-x1^2-x2^2)
     H <- rbind(c(2*x1*(2*x1^2-3)*E, 2*x2*(2*x1^2-1)*E),
                c(2*x2*(2*x1^2-1)*E, 2*x1*(2*x2^2-1)*E))
     
     ## sign of the sum of the eigenvalues of the Hessian
     z[i] <- sign(sum(eigen(H, symmetric=TRUE, only.values=TRUE)$values))
   }
   
   ## scale to {1,2}
   z[z == -1] <- 0
   C <- z + 1

   ## three classes or two
   if(threed) C[C==1 & X[,1] > 0] <- 3
   
   return(C)
 }

