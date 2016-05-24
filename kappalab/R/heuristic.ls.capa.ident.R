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

## Heuristic least squares capacity identification

##############################################################################

heuristic.ls.capa.ident <- function(n, mu, C, g, Integral="Choquet",maxiter = 500,
                                    alpha = 0.01, epsilon = 1e-6) {
    
    ## check n
    if (!(as.integer(n) == n))
        stop("wrong arguments")
    
    ## check mu
    if (!("capacity" %in% is(mu)))
        stop("wrong object mu")
    
    ## check C
    if (!(is.matrix(C) && dim(C)[2] == n)) 
        stop("wrong criteria matrix")
    
    ## check g
    if (!(is.vector(g) && length(g) == dim(C)[1])) 
        stop("wrong vector of global scores")

    ## check Integral
    if (!(Integral %in% c("Choquet","Sipos")))
        stop("wrong arguments")

    ## check alpha
    if (!(is.positive(alpha) && alpha <= 1))
        stop("wrong alpha value")

    ## check epsilon
    if (!(is.positive(epsilon) && epsilon <= 1e-2))
        stop("wrong epsilon value")
    
    ## number of alternatives
    n.alt <- dim(C)[1]

    ## heuristic least squares
    obj <-  .C("hlms", 
               as.integer(n),
               as.integer(Integral == "Choquet"),
               n.iter = as.integer(maxiter),
               mu = as.double(mu@data),
               as.integer(n.alt),
               as.double(t(C)),
               as.double(g),
               as.double(alpha),
               as.double(epsilon),
               error = double(1),
               PACKAGE="kappalab")

    ## solution
    mu <- as.game(mu)
    mu@data <- obj$mu

    if (is.monotone(mu))
        mu <- as.capacity(mu)
    else
        warning("The solution is not monotone")

    ## residuals
    residuals <- numeric(n.alt)
    for (i in 1:n.alt)
        if (Integral == "Choquet")
            residuals[i] <- g[i] - Choquet.integral(mu,C[i,])
        else
            residuals[i] <- g[i] - Sipos.integral(mu,C[i,])
    
    ## solution
    return(list(solution = mu, n.iter = obj$n.iter, residuals = residuals, mse = obj$error))
} 

##############################################################################
