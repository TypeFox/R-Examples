###############################################################################
## Classical optimal IC (optimal in sense of the Cramer-Rao bound)
###############################################################################
setMethod("optIC", signature(model = "L2ParamFamily", risk = "asCov"),
    function(model, risk){
        Curve <- as((trafo(model@param) %*% solve(model@FisherInfo)) %*% model@L2deriv, "EuclRandVariable")
        asCov <- trafo(model@param) %*% solve(model@FisherInfo) %*% t(trafo(model@param))

        modifyIC <- function(L2Fam, IC){ optIC(L2Fam, asCov()) }
        L2call <- model@fam.call
        L2call$trafo <- trafo(model)
        return(IC(
            name = paste("Classical optimal influence curve for", model@name), 
            CallL2Fam = L2call,
            Curve = EuclRandVarList(Curve),
            modifyIC = modifyIC,
            Risks = list(asCov = asCov, trAsCov = sum(diag(asCov))),
            Infos = matrix(c("optIC", "optimal IC in sense of Cramer-Rao bound"), 
                        ncol = 2, dimnames = list(character(0), c("method", "message")))))
    })
