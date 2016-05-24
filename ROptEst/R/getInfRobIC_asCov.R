###############################################################################
## get classical optimal IC
###############################################################################
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asCov", 
                                   neighbor = "ContNeighborhood"),
    function(L2deriv, risk, neighbor, Finfo, trafo, verbose = NULL){

            if(missing(verbose)|| is.null(verbose))
               verbose <- getRobAStBaseOption("all.verbose")

            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            
            b <- abs(as.vector(A))*max(abs(q(L2deriv)(1)),abs(q(L2deriv)(0)))
            
            asCov <- A %*% t(trafo)
            r <- neighbor@radius
            Risk <- list(asCov = asCov, 
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = asCov, normtype = NormType()),
                         asMSE = list(value = asCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))
            w <- new("HampelWeight")
            clip(w) <- b
            cent(w) <- 0
            stand(w) <- A
            weight(w) <- getweight(w, neighbor = neighbor,
                                   biastype = symmetricBias(),
                                   normW = NormType())

            return(list(A = A, a = 0, b = b, d = NULL, w = w, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "UnivariateDistribution", 
                                   risk = "asCov", 
                                   neighbor = "TotalVarNeighborhood"),
    function(L2deriv, risk, neighbor, Finfo, trafo, verbose = NULL){

            if(missing(verbose)|| is.null(verbose))
               verbose <- getRobAStBaseOption("all.verbose")

            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            b <- abs(as.vector(A))*(q(L2deriv)(1)-q(L2deriv)(0))
            a <- -abs(as.vector(A))*q(L2deriv)(0)
            asCov <- A %*% t(trafo)
            r <- neighbor@radius
            Risk <- list(asCov = asCov, 
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = NormType(), 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = asCov, normtype = NormType()),
                         asMSE = list(value = asCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))

            w <- new("BdStWeight")
            clip(w) <- c(0,b)+a
            stand(w) <- A
            weight(w) <- getweight(w, neighbor = neighbor,
                                   biastype = symmetricBias(),
                                   normW = NormType())

            return(list(A = A, a = -b/2, b = b, d = NULL, w = w, risk = Risk, info = info))
    })
setMethod("getInfRobIC", signature(L2deriv = "RealRandVariable", 
                                   risk = "asCov", 
                                   neighbor = "UncondNeighborhood"),
    function(L2deriv, risk, neighbor, Distr, Finfo, trafo, 
             QuadForm = diag(nrow(trafo)), verbose = NULL){

            if(missing(verbose)|| is.null(verbose))
               verbose <- getRobAStBaseOption("all.verbose")

            Cont <- is(neighbor,"ContNeighborhood")
            p <- nrow(trafo)
            if(! Cont && p>1)
                 stop("Not yet implemented")
            info <- c("optimal IC in sense of Cramer-Rao bound")
            A <- trafo %*% solve(Finfo)
            IC <- A %*% L2deriv
            if(is(Distr, "UnivariateDistribution")){
                lower <- ifelse(is.finite(q(Distr)(0)), q(Distr)(1e-8), q(Distr)(0))
                upper <- ifelse(is.finite(q(Distr)(1)), q(Distr)(1-1e-8), q(Distr)(1))
                x <- seq(from = lower, to = upper, length = 1e5)
                x <- x[x!=0] # problems with NaN=log(0)!
                ICx <- evalRandVar(IC, as.matrix(x))
                if(Cont)
                   b <- sqrt(max(colSums(ICx^2, na.rm = TRUE)))
                else{
                   b <- max(ICx)-min(ICx)
                }
            }else{
                b <- Inf # not yet implemented
            }
            
            asCov <- A %*% t(trafo)
            trAsCov <- sum(diag(QuadForm%*%asCov))
            r <- neighbor@radius
            nt <- if(identical(QuadForm,diag(nrow(trafo)))) NormType() else 
                     QFNorm(QuadForm = PosSemDefSymmMatrix(QuadForm))
            Risk <- list(asCov = asCov,  
                         asBias = list(value = b, biastype = symmetricBias(), 
                                       normtype = nt, 
                                       neighbortype = class(neighbor)), 
                         trAsCov = list(value = trAsCov, normtype = nt),
                         asMSE = list(value = trAsCov + r^2*b^2, 
                                      r = r,
                                      at = neighbor))
            if(Cont){
               w <- new("HampelWeight")
               clip(w) <- b
               cent(w) <- 0
               stand(w) <- A
            }else{
               w <- new("BdStWeight")
               clip(w) <- c(0,b)+as.numeric(A%*%z)
               stand(w) <- A
            }
            weight(w) <- getweight(w, neighbor = neighbor, biastype = symmetricBias(),
                                   normW = NormType())
            return(list(A = A, a = numeric(nrow(trafo)), b = b, d = NULL, w = w, risk = Risk,
                        info = info))
    })
