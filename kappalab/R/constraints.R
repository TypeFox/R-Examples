##############################################################################
#
# Copyright © 2005 Michel Grabisch and Ivan Kojadinovic    
#
# Ivan.Kojadinovic@polytech.univ-nantes.fr
#
# This software is a package for the statistical system GNU R:
# http://www.r-project.org 
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
#
##############################################################################

## Internal functions for handling the constraints of the quadratic programs

##############################################################################

## c <= C(a) - C(b) 

Choquet.preorder.constraint <- function(n, k, subsets, a, b, c) {
    
    if (!(is.positive(a) && is.positive(b) && is.positive(c) 
          && length(a) == n && length(b) == n && length(c) == 1))
        stop("wrong Choquet preorder constraint matrix")	
    
    A <- .C("Choquet_preorder_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.double(a), 
            as.double(b), 
            A = double(length(subsets)-1), 
            PACKAGE="kappalab")$A

    return(list(A = A, b = c))	
}

## Sh(i) >= Sh(j) + c

Shapley.preorder.constraint <- function(n, k, subsets, i, j, c) {
    
    if (!(i %in% 1:n && j %in% 1:n && i != j && length(c) == 1))
        stop("wrong Shapley preorder constraint matrix")	

    A <- .C("Shapley_preorder_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.integer(i-1), 
            as.integer(j-1), 
            A = double(length(subsets)-1), 
            PACKAGE="kappalab")$A

    return(list(A = A, b = c, r = 1 - c))	
}

## a <= Sh(i) <= b

Shapley.interval.constraint <- function(n, k, subsets, i, a, b) {
    
    if (!(i %in% 1:n && length(a) == 1 && length(b) == 1 &&
          a <= b && a >= 0 && b <= 1))
        stop("wrong Shapley interval constraint matrix")	

    A <- .C("Shapley_interval_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.integer(i-1), 
            A = double(length(subsets)-1), 
            PACKAGE="kappalab")$A

    return(list(A = A, b = a, r = b - a))	
}

## I(i1i2) >= I(j1j2) + c

interaction.preorder.constraint <- function(n, k, subsets, i1, i2,
                                                     j1, j2, c) {
    
    if (!(i1 %in% 1:n && i2 %in% 1:n && j1 %in% 1:n && j2 %in% 1:n
          && i1 != i2 && (i1 != j1 || i1 != j2 || i2 != j1
          || i2 != j2) && j1 != j2 && length(c) == 1))
        stop("wrong interaction preorder constraint matrix")	

    A <- .C("interaction_preorder_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.integer(i1-1), 
            as.integer(i2-1),
            as.integer(j1-1), 
            as.integer(j2-1), 
            A = double(length(subsets)-1), 
            PACKAGE="kappalab")$A

    return(list(A = A, b = c, r = 1 - c))	
}

## a <= I(ij) <= b

interaction.interval.constraint <- function(n, k, subsets, i, j, a, b) {
    
    if (!(i %in% 1:n && j %in% 1:n && i != j && length(a) == 1
          && length(b) == 1 && a <= b && a >= -1 && b <= 1))
        stop("wrong interaction interval constraint matrix")	

    A <- .C("interaction_interval_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.integer(i-1), 
            as.integer(j-1),
            A = double(length(subsets)-1), 
            PACKAGE="kappalab")$A

    return(list(A = A, b = a, r = b - a))	
}

## The partition {A1,...,Ap} is a inter-additive

inter.additive.partition.constraint <- function(n, k, subsets,
                                                partition) {
    n.var <- length(subsets)-1

    C <- .C("inter_additive_constraint", 
            as.integer(n),	
            as.integer(k), 
            as.integer(subsets),
            as.integer(partition),
            as.integer(max(partition)),
            C = double(n.var), 
            PACKAGE="kappalab")$C

    n.con <- sum(C)
    
    A <- matrix(0,n.con,n.var)

    i <- 1
    for (j in 1:n.var) 
        if (C[j] == 1) {

            A[i,j] <- 1
            i <- i + 1
        }

    return(list(A = A, b = rep(0,n.con), r = rep(0,n.con)))
}
