###############################################################################
## get classical optimal IC
###############################################################################
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- max(abs(as.vector(A)))*max(q(ErrorL2deriv)(1),abs(q(ErrorL2deriv)(0)))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*max(abs(q(Regressor)(1)), abs(q(Regressor)(0)))
                
            Risk <- list(asCov = A %*% t(trafo), asBias = b)

            return(list(A = A, a = numeric(nrow(trafo)), b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "UnivariateDistribution",
                                          risk = "asCov", 
                                          neighbor = "TotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- abs(as.vector(A))*(q(ErrorL2deriv)(1) - q(ErrorL2deriv)(0))
            b <- b*(abs(q(Regressor)(1)) + abs(q(Regressor)(0)))
            Risk <- list(asCov = A %*% t(trafo), asBias = b)

            return(list(A = A, a = -b/2, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- max(abs(as.vector(A)))*max(q(ErrorL2deriv)(1),abs(q(ErrorL2deriv)(0)))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*max(abs(q(Regressor)(1)), abs(q(Regressor)(0)))
            b.fct <- function(x){ b }
            body(b.fct) <- substitute({ b }, list(b = b))
            bfun <- RealRandVariable(Map = list(b.fct), 
                                     Domain = EuclideanSpace(dimension = dimension(img(Regressor))))
            Risk <- list(asCov = A %*% t(trafo), asBias = b*E(Regressor, neighbor@radiusCurve))

            return(list(A = A, a = numeric(nrow(trafo)), b = bfun, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- abs(as.vector(A))*(q(ErrorL2deriv)(1) - q(ErrorL2deriv)(0))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*(abs(q(Regressor)(1)) + abs(q(Regressor)(0)))
            b.fct <- function(x){ b }
            body(b.fct) <- substitute({ b }, list(b = b))
            bfun <- RealRandVariable(Map = list(b.fct), 
                                     Domain = EuclideanSpace(dimension = dimension(img(Regressor))))
            a.fct <- function(x){ -b/2 }
            body(a.fct) <- substitute({ -b/2 }, list(b = b))
            afun <- RealRandVariable(Map = list(a.fct), 
                                     Domain = EuclideanSpace(dimension = dimension(img(Regressor))))
            Risk <- list(asCov = A %*% t(trafo), asBias = b*E(Regressor, neighbor@radiusCurve))

            return(list(A = A, a = afun, b = bfun, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- max(abs(as.vector(A)))*max(q(ErrorL2deriv)(1),abs(q(ErrorL2deriv)(0)))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*max(abs(q(Regressor)(1)), abs(q(Regressor)(0)))
            Risk <- list(asCov = A %*% t(trafo), asBias = b)
            a.fct <- function(x){numeric(k)}
            body(a.fct) <- substitute({numeric(k)}, list(k = nrow(trafo)))
            Dom <- EuclideanSpace(dimension = dimension(img(Regressor)) + 1)
            a <- EuclRandVarList(EuclRandVariable(Map = list(a.fct), Domain = Dom, 
                                              dimension = trunc(nrow(trafo))))

            return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "Av2CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- max(abs(as.vector(A)))*max(q(ErrorL2deriv)(1),abs(q(ErrorL2deriv)(0)))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*max(abs(q(Regressor)(1)), abs(q(Regressor)(0)))
            Risk <- list(asCov = A %*% t(trafo), asBias = b)

            return(list(A = 1, z = 0, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "UnivariateDistribution", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "Av1CondTotalVarNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorL2derivDistrSymm, 
             RegSymm, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- max(abs(as.vector(A)))*abs(q(ErrorL2deriv)(1) - q(ErrorL2deriv)(0))
            if(is(Regressor, "UnivariateDistribution"))
                b <- b*(q(Regressor)(1) - q(Regressor)(0))
            Risk <- list(asCov = A %*% t(trafo), asBias = b)
            a.fct <- function(x){-b/2}
            body(a.fct) <- substitute({-b/2}, list(b = b))
            a <- RealRandVariable(Map = list(a.fct), Domain = img(Regressor))

            return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "ContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)

            if(is(ErrorDistr, "UnivariateDistribution")){
                lower <- ifelse(is.finite(q(ErrorDistr)(0)), q(ErrorDistr)(1e-8), q(ErrorDistr)(0))
                upper <- ifelse(is.finite(q(ErrorDistr)(1)), q(ErrorDistr)(1-1e-8), q(ErrorDistr)(1))
                x <- seq(from = lower, to = upper, length = 1e4)
                x <- x[x!=0] # problems with NaN=log(0)!
                b <- evalRandVar(ErrorL2deriv, as.matrix(x))^2
                b <- sqrt(max(colSums(b, na.rm = TRUE)))
            }else{
                b <- Inf # not yet implemented
            }
            asCov <- A %*% t(trafo)
            Risk <- list(asCov = asCov, asBias = b, trAsCov = sum(diag(asCov)))

            return(list(A = A, a = numeric(nrow(trafo)), b = b, d = NULL, risk = Risk, info = info))
    })
setMethod("getInfRobRegTypeIC", signature(ErrorL2deriv = "RealRandVariable", 
                                          Regressor = "Distribution",
                                          risk = "asCov", 
                                          neighbor = "Av1CondContNeighborhood"),
    function(ErrorL2deriv, Regressor, risk, neighbor, ErrorDistr, Finfo, trafo){
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)

            if(is(ErrorDistr, "UnivariateDistribution")){
                lower <- ifelse(is.finite(q(ErrorDistr)(0)), q(ErrorDistr)(1e-8), q(ErrorDistr)(0))
                upper <- ifelse(is.finite(q(ErrorDistr)(1)), q(ErrorDistr)(1-1e-8), q(ErrorDistr)(1))
                x <- seq(from = lower, to = upper, length = 1e4)
                x <- x[x!=0] # problems with NaN=log(0)!
                b <- evalRandVar(ErrorL2deriv, as.matrix(x))^2
                b <- sqrt(max(colSums(b, na.rm = TRUE)))
            }else{
                b <- Inf # not yet implemented
            }
            asCov <- A %*% t(trafo)
            Risk <- list(asCov = asCov, asBias = b, trAsCov = sum(diag(asCov)))
            a.fct <- function(x){numeric(k)}
            body(a.fct) <- substitute({numeric(k)}, list(k = nrow(trafo)))
            Dom <- EuclideanSpace(dimension = dimension(img(Regressor)) + 1)
            a <- EuclRandVarList(EuclRandVariable(Map = list(a.fct), Domain = Dom, 
                                                  dimension = trunc(nrow(trafo))))
                                              
            return(list(A = A, a = a, b = b, d = NULL, risk = Risk, info = info))
    })
