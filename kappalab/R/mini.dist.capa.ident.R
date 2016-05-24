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

## Minimum distance capacity identification

##############################################################################

## Constructs a Mobius.capacity object by means of a quadratic program 

mini.dist.capa.ident <- function(a, k, distance = "Choquet.coefficients",
                                      A.Choquet.preorder = NULL,
                                      A.Shapley.preorder = NULL,
                                      A.Shapley.interval = NULL,
                                      A.interaction.preorder = NULL,
                                      A.interaction.interval = NULL,
                                      A.inter.additive.partition = NULL,
                                      epsilon = 1e-6) {

    
    ## check a
    if (!("Mobius.game" %in% is(a)))
        stop("Object a is not of class Mobius.game")
    
    ## number of elements (criteria)
    n <- a@n

    ## check k 
    if (!(k %in% 1:n))
        stop("wrong arguments")

    # check distance
    if (!(distance %in% c("Choquet.coefficients","binary.alternatives","global.scores")))
        stop("wrong distance type")

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

    ## number of variables
    n.var <- binom.sum(n,k) - 1

    ## size of the vector representing a
    s.a <- binom.sum(n,a@k) - 1

    ## 2^n - 1
    pow.n <- 2^n - 1
    
    ## number of monotonicity constraints
    n.con <- n*2^(n-1)

    ## k power set or a@k power set in natural order
    power.set <-  .C("k_power_set", 
                   as.integer(n),
                   as.integer(max(k,a@k)),
                   subsets = integer(max(n.var+1,s.a+1)),
                   PACKAGE="kappalab")$subsets

    subsets <- power.set[1:(n.var+1)]
    subsets.a <- power.set[1:(s.a+1)]
    
    ## monotonicity constraints
    M <- .C("monotonicity_constraints", 
            as.integer(n), 
            as.integer(k),
            as.integer(subsets),
            M = integer(n.var * n.con),
            PACKAGE="kappalab")$M


    M <- matrix(M,n.con,n.var,byrow=TRUE)

    ## objective function
    if (distance == "Choquet.coefficients") {
        
        cat("Distance used: Choquet.coefficients \n\n")
        
        D.Shapley <- .C("objective_function_Choquet_coefficients", 
                        as.integer(n), 
                        D = double(n.con),
                        PACKAGE="kappalab")$D		
        
        Dmat <- t(M) %*% diag(D.Shapley) %*% M
        
        ## forming the second part of the objective function (dvec)
        M.a <- .C("monotonicity_constraints", 
                  as.integer(n), 
                  as.integer(a@k),
                  as.integer(subsets.a),
                   M = integer(s.a * n.con),
                  PACKAGE="kappalab")$M
        
        M.a <- matrix(M.a,n.con,s.a,byrow=TRUE)
        

        ## linear part of the objective function
        dvec <- a@data[-1] %*% t(M.a) %*% diag(D.Shapley) %*% M
        
    } else if (distance == "binary.alternatives") {

        cat("Distance used: binary.alternatives \n\n")
        
        B <- .C("objective_function_binary_alternatives", 
                as.integer(n), 
                as.integer(k),
                as.integer(subsets),
                B = integer(n.var * pow.n),
                PACKAGE="kappalab")$B
        
        B <- matrix(B,pow.n,n.var,byrow=TRUE)

        Dmat <- t(B) %*% B
        
        dvec <- zeta(a)@data[-1] %*% B 
        
    } else { # distance = "global.scores"

        cat("Distance used: global.scores \n\n")
        
        Q <- .C("objective_function_global_scores", 
                as.integer(n), 
                as.integer(max(a@k,k)),
                as.integer(k),
                as.integer(power.set),
                Q = double(max(s.a,n.var) * n.var),
                PACKAGE="kappalab")$Q
        
        Q <- matrix(Q,max(s.a,n.var),n.var,byrow=TRUE)
        
        Dmat <- Q[1:n.var,]

        dvec <- a@data[-1] %*% Q[1:s.a,]
        
    }

    ## the constraint matrix
    A <- M
    
    ## add the normalization constraint sum a(T) = 1
    A <- rbind(rep(1,n.var),A)
    bvec <- c(1,rep(epsilon,n.con))
    meq <- 1
    
    ## add the constraints relative to the preorder of the alternatives
    if (!is.null(A.Choquet.preorder)) {
        
        for (i in 1:dim(A.Choquet.preorder)[1]) {
            
            cpc <- Choquet.preorder.constraint(n,k,subsets,
                                               A.Choquet.preorder[i,][1:n],
                                               A.Choquet.preorder[i,][(n+1):(2*n)],
                                               A.Choquet.preorder[i,2*n+1])

            A <- rbind(A,cpc$A)
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
            
            A <- rbind(A,spc$A)
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
            A <- rbind(A,sic$A)
            bvec <- c(bvec,sic$b)
            ## - Sh(i) >= -b
            A <- rbind(A,-sic$A)
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
            
            A <- rbind(A,ipc$A)
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
            A <- rbind(A,iic$A)
            bvec <- c(bvec,iic$b)
            ## I(ij) <= b
            A <- rbind(A,-iic$A)
            bvec <- c(bvec,-(iic$b+iic$r))
         }	
    }

    ## add the constraints relative to the inter-addtive partition
    if (!is.null(A.inter.additive.partition)) {
        
            
        iapc <- inter.additive.partition.constraint(n,k,subsets,
                                                 A.inter.additive.partition)

        A <- rbind(iapc$A,A)
        bvec <- c(iapc$b,bvec)
        meq <- meq + length(iapc$b)
    }
    
    ## quadprog
    qp <- solve.QP(Dmat, dvec , t(A), bvec, meq = meq)

    return(list(solution = Mobius.capacity(c(0,qp$solution),n,k),
                value = qp$value, iterations = qp$iterations,
                iact = qp$iact))    
    
} 

##############################################################################
