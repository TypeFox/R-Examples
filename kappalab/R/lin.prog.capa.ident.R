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

## Capacity identification using the linear programming based approach
## of Marichal and Roubens

##############################################################################

## Constructs a Mobius.capacity object by means of a linear program 

lin.prog.capa.ident <- function(n, k,
                                      A.Choquet.preorder = NULL,
                                      A.Shapley.preorder = NULL,
                                      A.Shapley.interval = NULL,
                                      A.interaction.preorder = NULL,
                                      A.interaction.interval = NULL,
                                      A.inter.additive.partition = NULL,
                                      epsilon = 1e-6) {

    ## check n and k
    if (!(as.integer(n) == n && k %in% 1:n))
        stop("wrong arguments")

    ## check A.Choquet.preorder
    if (!((is.matrix(A.Choquet.preorder)
           && dim(A.Choquet.preorder)[2] == 2*n+1) 
          || is.null(A.Choquet.preorder)))
        stop("wrong Choquet preorder constraint matrix")
    
    ## check A.Shapley.preorder
    if (!((is.matrix(A.Shapley.preorder) && dim(A.Shapley.preorder)[2] == 3) 
          || is.null(A.Shapley.preorder)))
        stop("wrong Shapley preorder constraint matrix")

    ## check A.Shapley.interval
    if (!((is.matrix(A.Shapley.interval) && dim(A.Shapley.interval)[2] == 3) 
          || is.null(A.Shapley.interval)))
        stop("wrong Shapley interval constraint matrix")

    ## check A.interaction.preorder
    if (!((is.matrix(A.interaction.preorder)
           && dim(A.interaction.preorder)[2] == 5)
          || is.null(A.interaction.preorder)))
        stop("wrong interaction preorder constraint matrix")

    ## check A.interaction.interval
    if (!((is.matrix(A.interaction.interval)
           && dim(A.interaction.interval)[2] == 4) 
          || is.null(A.interaction.interval)))
        stop("wrong interaction interval constraint matrix")

    ## check A.inter.additive.partition
    if (!((is.numeric(A.inter.additive.partition)
           && sum(levels(factor(A.inter.additive.partition))
                  == 1:max(A.inter.additive.partition))
           == max(A.inter.additive.partition))
          || is.null(A.inter.additive.partition)))
        stop("wrong inter-additive partition")

    ## check epsilon
    if (!(is.positive(epsilon) && epsilon <= 1e-3))
        stop("wrong epsilon value")

    ## number of variables without the slack variable z/Epsilon
    n.var <- binom.sum(n,k) - 1

    ## number of monotonicity constraints
    n.con <- n*2^(n-1)

    ## k power set in natural order
    subsets <-  .C("k_power_set", 
                   as.integer(n),
                   as.integer(k),
                   subsets = integer(n.var+1),
                   PACKAGE="kappalab")$subsets

    ## monotonicity constraints
    A <- .C("monotonicity_constraints", 
            as.integer(n), 
            as.integer(k),
            as.integer(subsets),
            A = integer(n.var * n.con),
            PACKAGE="kappalab")$A


    A <- matrix(A,n.con,n.var,byrow=TRUE)
    A <- cbind(A,rep(0,n.con))
    ineqvec <- rep(">=",n.con)

    ## add the normalization constraint sum a(T) = 1
    A <- rbind(c(rep(1,n.var),0),A)
    ineqvec <- c("==",ineqvec)
    bvec <- c(1,rep(epsilon,n.con))
    
    ## add the constraints relative to the preorder of the alternatives
    if (!is.null(A.Choquet.preorder)) {
        
        for (i in 1:dim(A.Choquet.preorder)[1]) {
            
            cpc <- Choquet.preorder.constraint(n,k,subsets,
                                               A.Choquet.preorder[i,][1:n],
                                               A.Choquet.preorder[i,][(n+1):(2*n)],
                                               A.Choquet.preorder[i,2*n+1])

            A <- rbind(A,c(cpc$A,-1))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,cpc$b)	
        }	
    }
    
    ## add the constraints relative to the preorder of the criteria
    if (!is.null(A.Shapley.preorder)) {
        
        for (i in 1:dim(A.Shapley.preorder)[1]) {
            
            spc <- Shapley.preorder.constraint(n,k,subsets,
                                             A.Shapley.preorder[i,1],
                                             A.Shapley.preorder[i,2],
                                             A.Shapley.preorder[i,3])
            
            A <- rbind(A,c(spc$A,-1))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,spc$b)	
        }	
    }

    ## add the constraints relative to the importance of the criteria
    if (!is.null(A.Shapley.interval)) {
        
        for (i in 1:dim(A.Shapley.interval)[1]) {
            
            sic <- Shapley.interval.constraint(n,k,subsets,
                                               A.Shapley.interval[i,1],
                                               A.Shapley.interval[i,2],
                                               A.Shapley.interval[i,3])

            ## Sh(i) >= a
            A <- rbind(A,c(sic$A,0))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,sic$b)
            ## - Sh(i) >= -b
            A <- rbind(A,c(-sic$A,0))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,-(sic$b + sic$r))
            
        }	
    }

    ## add the constraints relative to the preorder of the interactions
    if (!is.null(A.interaction.preorder)) {
        
        for (i in 1:dim(A.interaction.preorder)[1]) {
            
            ipc <- interaction.preorder.constraint(n,k,subsets,
                                             A.interaction.preorder[i,1],
                                             A.interaction.preorder[i,2],
                                             A.interaction.preorder[i,3],
                                             A.interaction.preorder[i,4],
                                             A.interaction.preorder[i,5])
            
            A <- rbind(A,c(ipc$A,-1))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,ipc$b)	
        }	
    }

    ## add the constraints relative to the magnitude of the interaction
    if (!is.null(A.interaction.interval)) {
        
        for (i in 1:dim(A.interaction.interval)[1]) {
            
            iic <- interaction.interval.constraint(n,k,subsets,
                                                   A.interaction.interval[i,1],
                                                   A.interaction.interval[i,2],
                                                   A.interaction.interval[i,3],
                                                   A.interaction.interval[i,4])

            ## I(ij) >= a
            A <- rbind(A,c(iic$A,0))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,iic$b)
            ## I(ij) <= b
            A <- rbind(A,c(-iic$A,0))
            ineqvec <- c(ineqvec,">=")
            bvec <- c(bvec,-(iic$b+iic$r))
        }	
    }

    ## add the constraints relative to the inter-addtive partition
    if (!is.null(A.inter.additive.partition)) {
        
            
        iapc <- inter.additive.partition.constraint(n,k,subsets,
                                                 A.inter.additive.partition)

        A <- rbind(A,cbind(iapc$A,0))
        bvec <- c(iapc$b,bvec)
        ineqvec <- c(ineqvec,rep("==",length(iapc$b)))
    }
    

    ## lower bound of the Mobius transform
    lower.bound <- .C("Mobius_lower_bound", 
                      as.integer(n), 
                      as.integer(k),
                      as.integer(subsets),
                      lb = double(n.var),
                      PACKAGE="kappalab")$lb

    # change of variable so that the LP be in standard form
    bvec <- bvec - A %*% c(lower.bound,0)
    
    ## lpSolve
    lp.res <- lp("max",c(rep(0,n.var),1),A,ineqvec,bvec)
    print(lp.res)

    return(list(solution = Mobius.capacity(c(0,lp.res$solution[1:n.var]+lower.bound),n,k),
                value = lp.res$objval,lp.object = lp.res))    
    
} 

##############################################################################
