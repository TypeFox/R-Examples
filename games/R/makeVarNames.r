##
## INPUT:
## varNames: list of character vectors, each containing the variable names for
## one utility or variance equation
## prefixes: character vector containing names of utility equations (but not
## variance terms)
## utils: integer containing indices of utility equations (or any others where
## "equation:(Intercept)" should not be truncated to "equation")
## link: "logit" or "probit"
## sdterms: number of variance equations
##
## RETURN:
## varNames: character vector of variable names
## hasColon: logical vector indicating which utility/variance equations are not
## fixed to 0 or contain only a constant
## 
makeVarNames <- function(varNames, prefixes, utils, link, sdterms)
{
    vname <- if (link == "logit") "log(lambda" else "log(sigma"
    if (sdterms == 1L) {
        prefixes <- c(prefixes, paste(vname, ")", sep = ""))
    } else if (sdterms > 1L) {
        prefixes <- c(prefixes, paste(vname, 1:sdterms, ")", sep = ""))
    }

    ## no colon for any equation with no terms (i.e. fixed to 0), but exclude
    ## colon for intercept-only equations that aren't utility equations
    ## (e.g. scale parameters in the ultimatum model)
    utils <- 1:length(varNames) %in% utils
    noTerms <- sapply(varNames, function(x) length(x) == 0)
    onlyConstant <- sapply(varNames, function(x) all(x == "(Intercept)"))
    hasColon <- rep(TRUE, length(varNames))
    hasColon[noTerms] <- FALSE
    hasColon[onlyConstant & !utils] <- FALSE
    names(hasColon) <- prefixes

    for (i in seq_along(varNames)) {
        if (hasColon[i]) {
            varNames[[i]] <- paste(prefixes[i], varNames[[i]], sep = ":")
        } else {
            varNames[[i]] <- prefixes[i][length(varNames[[i]])]
        }
    }
    varNames <- unlist(varNames)
    ans <- list(varNames = varNames, hasColon = hasColon)
    return(ans)
}
