######################################################################
############################ Class sClere ############################
############################## Creation ##############################
######################################################################


### Class definition ###
setClass(
    Class = "sClere", 
    representation = representation(
        analysis = "character", 
        g = "numeric", 
        nbVarGroups = "integer", 
        algorithm = "character", 
        intercept = "numeric", 
        b = "numeric", 
        pi = "numeric", 
        sigma2 = "numeric", 
        gamma2 = "numeric", 
        likelihood = "numeric", 
        entropy = "numeric", 
        AIC = "numeric", 
        BIC = "numeric", 
        ICL = "numeric"
    ), 
    prototype = prototype(
        analysis = character(), 
        g = numeric(), 
        nbVarGroups = integer(), 
        algorithm = character(), 
        intercept = numeric(), 
        b = numeric(), 
        pi = numeric(), 
        sigma2 = numeric(), 
        gamma2 = numeric(), 
        likelihood = numeric(), 
        entropy = numeric(), 
        AIC = numeric(), 
        BIC = numeric(), 
        ICL = numeric()
    )
)


### Show ###
setMethod(f = "show", signature = "sClere", definition = function(object){
    nd = 4
    sep = "\t"
    cat("\t-------------------------------\n")
    cat("\t| CLERE | Yengo et al. (2013) |\n")
    cat("\t-------------------------------\n\n")
    cat("\tModel object for")
    cat(paste(" ", object@nbVarGroups, "groups of variables"))
    select <- switch(EXPR = object@analysis, 
        "fit" = {" ( user-specified )\n"}, 
        "bic" = {paste(" ( Selected using BIC criterion )\n")}, 
        "aic" = {paste(" ( Selected using AIC criterion )\n")}, 
        "icl" = {paste(" ( Selected using ICL criterion )\n")}
    )
    cat(select)
    cat("\n\t---\n")
    cat(paste("\tEstimated parameters using", object@algorithm, "algorithm are\n"))
    cat("\tintercept = ");cat(format(object@intercept, digits = nd), sep = sep);cat("\n")
    cat("\tb         = ");cat(format(object@b, digits = nd), sep = sep);cat("\n")
    cat("\tpi        = ");cat(format(object@pi, digits = nd), sep = sep);cat("\n")
    cat("\tsigma2    = ");cat(format(object@sigma2, digits = nd), sep = sep);cat("\n")
    cat("\tgamma2    = ");cat(format(object@gamma2, digits = nd), sep = sep);cat("\n")
    cat("\n\t---\n")
    cat(paste("\tLog-likelihood = ", format(object@likelihood, digits = nd), "\n"))
    cat(paste("\tEntropy        = ", format(object@entropy, digits = nd+1), "\n"))
    cat(paste("\tAIC            = ", format(object@AIC, digits = nd+1), "\n"))
    cat(paste("\tBIC            = ", format(object@BIC, digits = nd+1), "\n"))
    cat(paste("\tICL            = ", format(object@ICL, digits = nd+1), "\n\n"))
    return(invisible(object))
})


### Getteur ###
setMethod(f = "[", signature = "sClere", definition = function(x, i, j, drop){
    switch(EXPR = i, 
        "analysis" = {
            if (missing(j)) {
                return(x@analysis)
            } else {
                if (j>length(x@analysis)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@analysis[j])
                }
            }
        }, 
        "g" = {
            if (missing(j)) {
                return(x@g)
            } else {
                if (j>length(x@g)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@g[j])
                }
            }
        }, 
        "nbVarGroups" = {
            if (missing(j)) {
                return(x@nbVarGroups)
            } else {
                if (j>length(x@nbVarGroups)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@nbVarGroups[j])
                }
            }
        }, 
        "algorithm" = {
            if (missing(j)) {
                return(x@algorithm)
            } else {
                if (j>length(x@algorithm)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@algorithm[j])
                }
            }
        }, 
        "intercept" = {
            if (missing(j)) {
                return(x@intercept)
            } else {
                if (j>length(x@intercept)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@intercept[j])
                }
            }
        }, 
        "b" = {
            if (missing(j)) {
                return(x@b)
            } else {
                if (j>length(x@b)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@b[j])
                }
            }
        }, 
        "pi" = {
            if (missing(j)) {
                return(x@pi)
            } else {
                if (j>length(x@pi)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@pi[j])
                }
            }
        }, 
        "sigma2" = {
            if (missing(j)) {
                return(x@sigma2)
            } else {
                if (j>length(x@sigma2)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@sigma2[j])
                }
            }
        }, 
        "gamma2" = {
            if (missing(j)) {
                return(x@gamma2)
            } else {
                if (j>length(x@gamma2)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@gamma2[j])
                }
            }
        }, 
        "likelihood" = {
            if (missing(j)) {
                return(x@likelihood)
            } else {
                if (j>length(x@likelihood)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@likelihood[j])
                }
            }
        }, 
        "entropy" = {
            if (missing(j)) {
                return(x@entropy)
            } else {
                if (j>length(x@entropy)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@entropy[j])
                }
            }
        }, 
        "AIC" = {
            if (missing(j)) {
                return(x@AIC)
            } else {
                if (j>length(x@AIC)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@AIC[j])
                }
            }
        }, 
        "BIC" = {
            if (missing(j)) {
                return(x@BIC)
            } else {
                if (j>length(x@BIC)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@BIC[j])
                }
            }
        }, 
        "ICL" = {
            if (missing(j)) {
                return(x@ICL)
            } else {
                if (j>length(x@ICL)) {
                    stop("[sClere:get] indice out of limits")
                } else {
                    return(x@ICL[j])
                }
            }
        }, 
        stop("[sClere:get] ", i, " is not a \"sClere\" slot")
    )
})


### Setteur ###
setMethod(f = "[<-", signature = "sClere", definition = function(x, i, j, value){
    switch(EXPR = i, 
         "analysis" = {
            if (missing(j)) {
                x@analysis <- value
            } else {
                if (j>length(x@analysis)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@analysis[j] <- value
                }
            }
        }, 
         "g" = {
            if (missing(j)) {
                x@g <- value
            } else {
                if (j>length(x@g)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@g[j] <- value
                }
            }
        }, 
         "nbVarGroups" = {
            if (missing(j)) {
                x@nbVarGroups <- value
            } else {
                if (j>length(x@nbVarGroups)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@nbVarGroups[j] <- value
                }
            }
        }, 
         "algorithm" = {
            if (missing(j)) {
                x@algorithm <- value
            } else {
                if (j>length(x@algorithm)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@algorithm[j] <- value
                }
            }
        }, 
         "intercept" = {
            if (missing(j)) {
                x@intercept <- value
            } else {
                if (j>length(x@intercept)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@intercept[j] <- value
                }
            }
        }, 
         "b" = {
            if (missing(j)) {
                x@b <- value
            } else {
                if (j>length(x@b)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@b[j] <- value
                }
            }
        }, 
         "pi" = {
            if (missing(j)) {
                x@pi <- value
            } else {
                if (j>length(x@pi)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@pi[j] <- value
                }
            }
        }, 
         "sigma2" = {
            if (missing(j)) {
                x@sigma2 <- value
            } else {
                if (j>length(x@sigma2)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@sigma2[j] <- value
                }
            }
        }, 
         "gamma2" = {
            if (missing(j)) {
                x@gamma2 <- value
            } else {
                if (j>length(x@gamma2)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@gamma2[j] <- value
                }
            }
        }, 
         "likelihood" = {
            if (missing(j)) {
                x@likelihood <- value
            } else {
                if (j>length(x@likelihood)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@likelihood[j] <- value
                }
            }
        }, 
         "entropy" = {
            if (missing(j)) {
                x@entropy <- value
            } else {
                if (j>length(x@entropy)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@entropy[j] <- value
                }
            }
        }, 
         "AIC" = {
            if (missing(j)) {
                x@AIC <- value
            } else {
                if (j>length(x@AIC)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@AIC[j] <- value
                }
            }
        }, 
         "BIC" = {
            if (missing(j)) {
                x@BIC <- value
            } else {
                if (j>length(x@BIC)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@BIC[j] <- value
                }
            }
        }, 
         "ICL" = {
            if (missing(j)) {
                x@ICL <- value
            } else {
                if (j>length(x@ICL)) {
                    stop("[sClere:set] indice out of limits")
                } else {
                    x@ICL[j] <- value
                }
            }
        }, 
        stop("[sClere:set] ", i, " is not a \"sClere\" slot")
    )
    validObject(x)
    return(invisible(x))
})

