getL2normL2deriv <-
        function(aFinfo, cent, ...){sqrt(aFinfo+cent^2)}

setMethod("getL1normL2deriv", signature(L2deriv = "UnivariateDistribution"),
    function(L2deriv, cent, ...){
        return(-2*m1df(L2deriv, cent) +cent*(2*p(L2deriv)(cent)-1))
    })

setMethod("getL1normL2deriv", signature(L2deriv = "RealRandVariable"),
    function(L2deriv, cent, stand, Distr, normtype, ...){
        integrandG <- function(x, L2, stand, cent){
            X <- evalRandVar(L2, as.matrix(x))[,,1] - cent
            Y <- apply(X, 2, "%*%", t(stand))
            res <- fct(normtype)(Y)
            return((res > 0)*res)
        }

        return(E(object = Distr, fun = integrandG, L2 = L2deriv,
                  stand = stand, cent = cent, useApply = FALSE))
    })
