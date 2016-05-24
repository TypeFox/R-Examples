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
################################################################################

## Unsupervised capacity identification from data based on entropy measures

##############################################################################

## Internal functions

## Converts a subset represented as an integer (binary notation) 
## to a numeric
binary2subset <- function(n,b) {
    
    result <- .C("binary2subsetR", 
                 as.integer(n),
                 as.integer(b),
                 subset = integer(n), 
                 len = integer(1),
                 PACKAGE="kappalab")
    
    result$subset[1:result$len]
}

shannon.function <- function(x) {
    if (x > 0)
        return(-x * log(x))
    else
        return(0)
}

shannon.entropy <- function(f,...) {
    
    sum(sapply(f,shannon.function))
}

renyi.entropy <- function(f,alpha) {
    
    1/(1-alpha)*log(sum(f^alpha))
}

havrda.charvat.entropy <- function(f,beta) {
    
    1/(1-beta)*(sum(f^beta)-1)
}

##############################################################################

##Constructs a capacity from discretized profiles
## and using parametric entropy measures  

entropy.capa.ident <- function(d,entropy = "renyi",parameter = 1)
{
    if (!is.data.frame(d))
        stop("wrong arguments")
    
    n <- length(d)

    for(i in 1:n)
        if (!is.factor(d[[i]]))
            stop("wrong arguments")

    if (!(entropy %in% c("renyi","havrda.charvat")))
        stop("wrong arguments")	

    if (!(as.double(parameter) && length(parameter) == 1))
        stop("wrong arguments")	

    if (parameter == 1)
        entropy.measure <- shannon.entropy
    else if (entropy == "renyi")
        entropy.measure <- renyi.entropy
    else 
        entropy.measure <- havrda.charvat.entropy

    ## multidimensional contingency table
    t <- table(d)
    
    ## frequency table
    f <- t/sum(t)

    ## unsupervised capacity identification

    mu <- numeric(2^n)
    
    for (i in 2:2^n)
        mu[i] <- entropy.measure(margin.table(f,binary2subset(n,i-1)),
                                 parameter)

    subsets <-  .C("k_power_set", 
                   as.integer(n),
                   as.integer(n),
                   subsets = integer(2^n),
                   PACKAGE="kappalab")$subsets

    new("capacity", data = mu/mu[2^n], subsets = subsets, n = n)
}

#############################################################################
