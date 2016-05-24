# Copyright (C) 2009 
# Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Kim-Anh Le Cao, French National Institute for Agricultural Research and 
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


vip <-
function(object)
{

    #-- initialisation des matrices --#
    W = object$loadings$X
    H = object$ncomp
    q = ncol(object$Y)
    p = ncol(object$X)
    VIP = matrix(0, nrow = p, ncol = H)
     
    cor2 = cor(object$Y, object$variates$X, use = "pairwise")^2
    cor2 = as.matrix(cor2, nrow = q)
     
    VIP[, 1] = W[, 1]^2
     
    if (H > 1) {
        for (h in 2:H) {
            if (q == 1) {
                Rd = cor2[, 1:h] 
                VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
            }
            else {
                Rd = apply(cor2[, 1:h], 2, sum)
                VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
            }
        }
    }
     
    #-- valeurs sortantes --#
    VIP = sqrt(p * VIP)
    rownames(VIP) = rownames(W)
    colnames(VIP)= paste("comp", 1:H)
     
    return(invisible(VIP))
}
