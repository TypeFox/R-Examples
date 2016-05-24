  ###################################################################
  #
  # This function is part of WACSgen V1.0
  # Copyright © 2013,2014,2015, D. Allard, BioSP,
  # and Ronan Trépos MIA-T, INRA
  #
  # This program is free software; you can redistribute it and/or
  # modify it under the terms of the GNU General Public License
  # as published by the Free Software Foundation; either version 2
  # of the License, or (at your option) any later version.
  #
  # This program is distributed in the hope that it will be useful,
  # but WITHOUT ANY WARRANTY; without even the implied warranty of
  # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  # GNU General Public License for more details. http://www.gnu.org 
  #
  ###################################################################
  sqrtm = function(M) {
  ###############################################
  #
  # Computes the square root of a p.d. matrix
  #
  # ARGUMENT 
  #     M: square matrix
  #
  # VALUE
  #    MS: square matrix, such that MS%*%MS = M
  #
  ###############################################
    if (nrow(M) != ncol(M)) stop("[sqrtm] the matrix should be a squared matrix")
    
    x1 = eigen(M)
    if (any(x1$values < 0)) {
      warning(paste("[sqrtm] the matrix should be positive definite:", min(x1$values)));
    }
    x1s  = sqrt(diag(x1$values))
    u1   = x1$vectors
    u1i  = solve(u1)
    MS   = u1 %*% x1s %*% u1i
    return(MS)
  }


  def.pos = function(M){
  ###############################################
  #
  # Find an approximation of M that is definite positive
  #
  # ARGUMENT 
  #     M   : square matrix
  #
  # VALUE
  #     Mdp : a square def.pos matrix
  #
  ###############################################
  if (nrow(M) != ncol(M)) stop("[def.pos] the matrix should be a squared matrix")
  
  N   = nrow(M)
  Mdp = M
  M.eigen = eigen(M)
  if ( M.eigen$values[N] < 0 ){
    dd = M.eigen$values
    neg=1
    for (i in 1:N) {
      if (dd[i] > 0) ddmin = dd[i]
      if (dd[i]  < 0){
        dd[i] = ddmin*(0.1^neg)
        neg=neg+1
      }
    } 
    Mdp = M.eigen$vectors%*%diag(dd)%*%t(M.eigen$vectors)
  }else{}
  return(Mdp)
  }
