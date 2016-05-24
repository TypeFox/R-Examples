getBoundedIC <- function(L2Fam, D=trafo(L2Fam@param)){
        FI <- FisherInfo(L2Fam)
        bm <- sum(diag(solve(FI)))
        w <- new("BoundedWeight", clip = bm, weight = function(x){
                   norm0 <- EuclideanNorm(as.matrix(x))
                   ind2 <- (norm0 < bm/2)
                   norm1 <- ind2*bm/2 + (1-ind2)*norm0
                   ind1 <- (norm0 < bm)
                   return(ind1 + (1-ind1)*bm/norm1)})

        dims <- length(L2Fam@param)

        L2deriv <- as(diag(dims) %*% L2Fam@L2deriv, "EuclRandVariable")

        ICfct <- vector(mode = "list", length = dims)
        L.fct <- function(x) evalRandVar(L2deriv,x)

        for(i in 1:dims){
                ICfct[[i]] <- function(x){}
                body(ICfct[[i]]) <- substitute({ Yi(x)*w(L(x)) },
                                                 list(Yi = L2deriv@Map[[i]],
                                                      L = L.fct,
                                                      w = weight(w)))
            }
        L2w <- EuclRandVariable(Map = ICfct, Domain = L2deriv@Domain,
                                         Range = L2deriv@Range)
        D1 <- L2Fam@distribution

        cent <- E(D1,L2w)
        L2w0 <- L2w - cent

        E1 <- matrix(E(D1, L2w0 %*% t(L2deriv-cent)), dims, dims)
        stand <- as.matrix(D %*% solve(E1, generalized = TRUE))
        return(as(stand %*% L2w0, "EuclRandVariable"))
        }
