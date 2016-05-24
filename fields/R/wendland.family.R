# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
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
Wendland2.2 <- function(d, theta = 1) {
    # Cari's test function  with explicit form  for d=2 k=2
    # taper range is 1.0
    d <- d/theta
    if (any(d < 0)) 
        stop("d must be nonnegative")
    return(((1 - d)^6 * (35 * d^2 + 18 * d + 3))/3 * (d < 1))
}
#
# the monster
#
"wendland.cov" <- function(x1, x2=NULL, theta = 1, V = NULL, 
    k = 2, C = NA, marginal = FALSE, Dist.args = list(method = "euclidean"), 
    spam.format = TRUE, derivative = 0,  verbose = FALSE) {
    #
    #   if marginal variance is needed
    #  this is a quick return
    if (marginal) {
        return(rep(1, nrow(x1)))
    }
    #  the rest of the possiblities require some computing
    # setup the two matrices of locations
    #
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
        }
    if( is.null( x2) ) {
    	 x2<- x1}  
    if (!is.matrix(x2) ) {
        x2 <- as.matrix(x2)
        }
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    # logical to figure out if this is great circle distance or not
    # great circle needs to handled specially due to how things are scaled.
    great.circle <- Dist.args$method == "greatcircle"
    # derivatives are tricky for great circle and other distances and have not been implemented ...
    if (Dist.args$method != "euclidean" & derivative > 0) {
        stop("derivatives not supported for this distance metric")
    }
    # catch bad theta format
    if (length(theta) > 1) {
        stop("theta as a matrix or vector has been depreciated")
    }
    # catch using V with great circle
    if (!is.null(V) & great.circle) {
        stop("V is not implemented with great circle distance")
    }
    if (!is.null(V)) {
        if (theta != 1) {
            stop("can't specify both theta and V!")
        }
        x1 <- x1 %*% t(solve(V))
        x2 <- x2 %*% t(solve(V))
    }
    # if great circle distance set the delta cutoff to be in scale of angular latitude.
    # also figure out if scale is in miles or kilometers
    if (great.circle) {
        miles <- ifelse(is.null(Dist.args$miles), TRUE, Dist.args$miles)
        delta <- (180/pi) * theta/ifelse(miles, 3963.34, 6378.388)
    }
    else {
        delta <- theta
    }
    if (verbose) {
        print(delta)
    }
    # once scaling is done taper is applied with default range of 1.0
    # find polynomial coeffients that define
    # wendland on [0,1]
    #  d dimension and  k is the order
    #  first  find sparse matrix of Euclidean distances
    #  ||x1-x2||**2 (or any other distance that may be specified by
    #   the method component in Dist.args
    #  any distance beyond delta is set to zero -- anticipating the
    # tapering to zero by the Wendland.
    #
    sM <- do.call("nearest.dist", c(list(x1, x2, delta = delta, 
        upper = NULL), Dist.args))
    # scale distances by theta
    # note: if V is passed then theta==1 and all the scaling should be done with the V matrix.
    # there are two possible actions listed below:
    # find  Wendland cross covariance matrix
    # return either in sparse or matrix format
    if (is.na(C[1])) {
        sM@entries <- Wendland(sM@entries/theta, k = k, dimension = d)
        if (!spam.format) {
            return(as.matrix(sM))
        }
        else {
            return(sM)
        }
    }
    else {
        #
        # multiply cross covariance matrix by the matrix C where
        # columns are usually the  'c' coefficients
        #  note multiply happens in spam format
        #
        if (derivative == 0) {
            sM@entries <- Wendland(sM@entries/theta, k = k, dimension = d)
            return(sM %*% C)
        }
        else {
            #        otherwise evaluate partial derivatives with respect to x1
            #        big mess of code and an explicit for loop!
            #         this only is for euclidean distance
            if (is.matrix(C)) {
                if (ncol(C) > 1) {
                  stop("C should be a vector")
                }
            }
            L <- length(coef)
            #         loop over dimensions and accumulate partial derivative matrix.
            tempD <- sM@entries
            tempW <- Wendland(tempD/theta, k = k, dimension = d, 
                derivative = derivative)
            # loop over dimensions and knock out each partial accumulate these in
            # in temp
            temp <- matrix(NA, ncol = d, nrow = n1)
            # Create rowindices vector
            sMrowindices <- rep(1:n1, diff(sM@rowpointers))
            for (kd in 1:d) {
                #
                #            Be careful if the distance (tempD) is close to zero.
                #            Note that the x1 and x2 are in transformed ( V inverse) scale
                #
                #
                sM@entries <- ifelse(tempD == 0, 0, (tempW * 
                  (x1[sMrowindices, kd] - x2[sM@colindices, kd])/(theta * 
                  tempD)))
                #
                # accumlate the new partial
                temp[, kd] <- sM %*% C
            }
            # transform back to original coordinates.
            if (!is.null(V)) {
                temp <- temp %*% t(solve(V))
            }
            return(temp)
        }
    }
    # should not get here!
}
#
#
#
Wendland2.2 <- function(d, theta = 1) {
    # Cari Kaufman's test case with explicit form  for d=2 k=2
    # taper range is 1.0
    d <- d/theta
    if (any(d < 0)) 
        stop("d must be nonnegative")
    return(((1 - d)^6 * (35 * d^2 + 18 * d + 3))/3 * (d < 1))
}
############## basic evaluation of Wendland and its derivatives.
###########################
# n: Wendland interpolation matrix is positive definite on R^n, i.e.  n is
# the dimension of the locations.
# k: Wendland function is 2k times continuously
# differentiable.
# The proofs can be found in the work of Wendland(1995).
#  H. Wendland. Piecewise polynomial , positive definite and compactly supported radial
#  functions of minimal degree. AICM 4(1995), pp 389-396.
#########################################
## top level function:
Wendland = function(d, theta = 1, dimension, k, derivative = 0, 
    phi = NA) {
    if (!is.na(phi)) {
        stop("phi argument has been depreciated")
    }
    if (any(d < 0)) {
        stop("d must be nonnegative")
    }
    # find scaling so that function at zero is 1.
    scale.constant <- wendland.eval(0, n = dimension, k, derivative = 0)
    # adjust by theta
    if (derivative > 0) {
        scale.constant <- scale.constant * (theta^(derivative))
    }
    # scale distances by theta.
    if( theta!=1){
         d <- d/theta}
    # at this point d the distances shouls be scaled so that
    # covariance is zero beyond 1
    if( (k==2)& (dimension==2) & (derivative==0)){
      ((1 - d)^6 * (35 * d^2 + 18 * d + 3))/3 * (d < 1)}
    else{
    ifelse(d < 1, wendland.eval(d, n = dimension, k, derivative)/scale.constant, 
        0)
    }
  }

####################
# [M] = wm(n, k)
# Compute the matrix coeficient in Wendland(1995)
# Input:
#\tn: Wendland interpolation matrix is positive definite on R^n
# \tk: Wendland function is 2k times continuously differentiable
####################
Wendland.beta = function(n, k) {
    l = floor(n/2) + k + 1
    M = matrix(0, nrow = k + 1, ncol = k + 1)
    #
    # top corner is 1
    #
    M[1, 1] = 1
    #
    # Compute across the columns and down the rows, filling out upper triangle of M (including diagonal). The indexing is done from 0, thus we have to adjust by +1 when accessing our matrix element.
    #
    if (k == 0) {
        stop
    }
    else {
        for (col in 0:(k - 1)) {
            #
            # Filling out the col+1 column
            #
            # As a special case, we need a different formula for the top row
            #
            row = 0
            beta = 0
            for (m in 0:col) {
                beta = beta + M[m + 1, col + 1] * fields.pochdown(m + 
                  1, m - row + 1)/fields.pochup(l + 2 * col - 
                  m + 1, m - row + 2)
            }
            M[row + 1, col + 2] = beta
            #
            # Now do the rest of rows
            #
            for (row in 1:(col + 1)) {
                beta = 0
                for (m in (row - 1):col) {
                  beta = beta + M[m + 1, col + 1] * fields.pochdown(m + 
                    1, m - row + 1)/fields.pochup(l + 2 * col - 
                    m + 1, m - row + 2)
                }
                M[row + 1, col + 2] = beta
            }
        }
    }
    M
}
########################################
# [phi] = wendland.eval(r, n, k, derivative).
# Compute the compacted support basis function in Wendland(1995).
# Input:
#\tr: a scalar representing the distance between locations. r should be scaled into [0,1] beforehand.
# \tn: Wendland interpolation matrix is positive definite on R^n. Or, we could say n is the dimension of the locations.
# \tk: Wendland function is 2k times continuously differentiable.
#\tderivative: the derivative of wendland function.
# Output:
#\tphi: a scalar evaluated by the Wendland function at distance r.
# example:
#\tr = 0.5
#\tphi = wendland.eval(r, 2, 1,derivative = 1 )
# The proofs can be found in the work of Wendland(1995).
# H. Wendlamd. Piecewise polynomial , positive definite and compactly supported radial functions of minimal degree. AICM 4(1995), pp 389-396.
#########################################
wendland.eval = function(r, n, k, derivative = 0) {
    #
    # check if the distances are between [0,1]
    #
    beta = Wendland.beta(n, k)
    l = floor(n/2) + k + 1
    if (derivative == 0) {
        #
        # first evaluate outside for loop with  m =0
        phi = beta[1, k + 1] * (1 - r)^(l + 2 * k)
        # now accumulate terms for other m values up to k
        for (m in 1:k) {
            phi = phi + beta[m + 1, k + 1] * r^m * (1 - r)^(l + 
                2 * k - m)
        }
    }
    else {
        # evaluate derivative note use of symbolic differtiation.
        f.my = expression((1 - r)^(l + 2 * k))
        f.deriv = fields.D(f.my, "r", order = derivative)
        f.eval = eval(f.deriv)
        phi = beta[1, k + 1] * f.eval
        for (m in 1:k) {
            f.my = expression(r^m * (1 - r)^(l + 2 * k - m))
            f.deriv = fields.D(f.my, "r", order = derivative)
            f.eval = eval(f.deriv)
            phi = phi + beta[m + 1, k + 1] * f.eval
        }
    }
    phi
}
#######################
# [n] = fields.pochup(q, k)
# Calculate the Pochhammer symbol for rising factorial q(q+1)(q+2)...(q+k-1)
#######################
fields.pochup = function(q, k) {
    n = q
    if (k == 0) {
        n = 1
    }
    else {
        for (j in 1:(k - 1)) {
            if ((k - 1) < 1) {
                stop
            }
            else {
                n = n * (q + j)
            }
        }
    }
    n
}
#########################
# [n] = fields.pochdown(q, k)
# Calculate the Pochhammer symbol for falling factorial q(q-1)(q-2)...(q-k+1)
#########################
fields.pochdown = function(q, k) {
    n = q
    if (k == 0) {
        n = 1
    }
    else {
        for (j in 1:(k - 1)) {
            if ((k - 1) < 1) {
                stop
            }
            else {
                n = n * (q - j)
            }
        }
    }
    n
}
#############################
# fields.D(f,name = x,order = n) forms the n-th derivative of function f with respect to the variable  x
################################
fields.D = function(f, name, order = 1) {
    if (order < 1) {
        stop("'order' must be >= 1")
    }
    if (order == 1) {
        d = D(f, name)
    }
    else {
        fields.D(D(f, name), name, order - 1)
    }
}
