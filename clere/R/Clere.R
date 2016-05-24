######################################################################
############################ Class Clere #############################
############################## Creation ##############################
######################################################################

### Class definition ###
setClass(
    Class = "Clere", 
    representation = representation(
        y = "numeric", 
        x = "matrix", 
        n = "integer", 
        p = "integer", 
        g = "numeric", 
        nItMC = "numeric", 
        nItEM = "numeric", 
        nBurn = "numeric", 
        dp = "numeric", 
        nsamp = "numeric", 
        sparse = "logical",  
        analysis = "character",  
        algorithm = "character", 
        initialized = "logical", 
        maxit = "numeric", 
        tol = "numeric", 
        seed = "integer", 
        b = "numeric", 
        pi = "numeric", 
        sigma2 = "numeric", 
        gamma2 = "numeric", 
        intercept = "numeric", 
        likelihood = "numeric", 
        entropy = "numeric", 
        P = "matrix", 
        theta = "matrix", 
        Zw = "matrix", 
        Bw = "matrix", 
        Z0 = "numeric", 
        message = "character"
    ), 
    prototype = prototype(
        y = numeric(), 
        x = matrix(, nrow = 0, ncol = 0), 
        n = integer(), 
        p = integer(), 
        g = numeric(), 
        nItMC = numeric(), 
        nItEM = numeric(), 
        nBurn = numeric(), 
        dp = numeric(), 
        nsamp = numeric(), 
        sparse = logical(),  
        analysis = character(),  
        algorithm = character(), 
        initialized = logical(), 
        maxit = numeric(), 
        tol = numeric(), 
        seed = integer(), 
        b = numeric(), 
        pi = numeric(), 
        sigma2 = numeric(), 
        gamma2 = numeric(), 
        intercept = numeric(), 
        likelihood = numeric(), 
        entropy = numeric(), 
        P = matrix(, nrow = 0, ncol = 0), 
        theta = matrix(, nrow = 0, ncol = 0), 
        Zw = matrix(, nrow = 0, ncol = 0), 
        Bw = matrix(, nrow = 0, ncol = 0), 
        Z0 = numeric(), 
        message = character()
    )
)


### Constructor ###
### Show ###
setMethod(f = "show", signature = "Clere", definition = function(object) {
  showSlot <- function(slot) {
    sNames <- gsub("^[^@]*@(.*)", "\\1", slot)
    eSlot <- eval(parse(text = slot))
    tmp <- switch(EXPR = class(eSlot),
                  "matrix" = {
                    cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
                    if (all(dim(eSlot)==0)) {
                      cat("NA")
                    } else {
                      cat("\n")
                      nrowShow <- seq(min(5, nrow(eSlot)))
                      ncolShow <- seq(min(5, ncol(eSlot)))
                      shortObject <- eSlot[nrowShow, ncolShow]
                      if (is.null(rownames(shortObject))) {
                        rownames(shortObject) <- seq(nrow(shortObject))
                      }
                      if (is.null(colnames(shortObject))) {
                        colnames(shortObject) <- seq(ncol(shortObject))
                      }
                      resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
                      if (nrow(shortObject)!=nrow(eSlot)) {
                        resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function(iCol) {paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")})))
                      }
                      if (ncol(shortObject)!=ncol(eSlot)) {
                        resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat)-1)))
                      }
                      cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
                    }
                    cat("\n")
                  },
                  "data.frame" = {
                    cat(paste0(" ~ ", sNames, " : [", nrow(eSlot), "x", ncol(eSlot), "]", collapse = ""))
                    if (all(dim(eSlot)==0)) {
                      cat("NA")
                    } else {
                      cat("\n")
                      nrowShow <- seq(min(5, nrow(eSlot)))
                      ncolShow <- seq(min(5, ncol(eSlot)))
                      shortObject <- eSlot[nrowShow, ncolShow]
                      if (is.null(rownames(shortObject))) {
                        rownames(shortObject) <- seq(nrow(shortObject))
                      } else {}
                      if (is.null(colnames(shortObject))) {
                        colnames(shortObject) <- seq(ncol(shortObject))
                      } else {}
                      resFormat <- format(cbind(c("", rownames(shortObject)), rbind(colnames(shortObject), format(shortObject, digits = 4))), justify = "centre")
                      if (nrow(shortObject)!=nrow(eSlot)) {
                        resFormat <- rbind(resFormat, c(".", sapply(seq(colnames(shortObject)), function(iCol) {paste0(rep(".", nchar(resFormat[1, 1])), collapse = "")})))
                      } else {}
                      if (ncol(shortObject)!=ncol(eSlot)) {
                        resFormat <- cbind(resFormat, c(".", rep(paste0(rep(".", nchar(resFormat[1, 1])), collapse = ""), nrow(resFormat)-1)))
                      } else {}
                      cat(paste0("     ", apply(format(resFormat, justify = "centre"), 1, paste, collapse = " "), "\n", collapse = ""))
                    }
                    cat("\n")
                  },
                  "numeric" = {
                    cat(paste0(" ~ ", sNames, " : ", collapse = ""))
                    if (length(eSlot) == 0) {
                      cat("NA")
                    } else {
                      if (length(eSlot)>1) {
                        cat(paste0("[", length(eSlot), "] ", paste0(format(head(eSlot), digits = 4), collapse = " ")))
                      } else {
                        cat(format(eSlot, digits = 4))
                      }
                    }
                    cat("\n")
                  },
                  "character" = {
                    cat(paste0(" ~ ", sNames, " : ", collapse = ""))
                    if (length(eSlot) == 0) {
                      cat("NA")
                    } else {
                      if (length(eSlot)>1) {
                        cat("[", length(eSlot), "] \"", paste0(head(eSlot), collapse = "\" \""), "\"", sep = "")
                      } else {
                        cat(paste0("\"", eSlot, "\""))
                      }
                    }
                    cat("\n")
                  },
                  {
                    cat(paste0(" ~ ", sNames, " : ", collapse = ""))
                    if (length(eSlot) == 0) {
                      cat("NA")
                    } else {
                      if (length(eSlot)>1) {
                        cat(paste0("[", length(eSlot), "] ", paste0(head(eSlot), collapse = " ")))
                      } else {
                        cat(eSlot)
                      }
                    }
                    cat("\n")
                  }
                  )
    return(invisible())
  }
  showObject <- function(object) {
    cat("	~~~ Class:", class(object), "~~~\n")
    sNames <- paste0("object@", slotNames(object))
    trash <- sapply(sNames, showSlot)
    return(invisible())
  }
  showObject(object)
  return(invisible(object))
})


### Getteur ###
setMethod(f = "[", signature = "Clere", definition = function(x, i, j, drop) {
    switch(EXPR = i, 
        "y" = {
            if (missing(j)) {
                return(x@y)
            } else {
                if (j>length(x@y)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@y[j])
                }
            }
        }, 
        "x" = {return(x@x)}, 
        "n" = {
            if (missing(j)) {
                return(x@n)
            } else {
                if (j>length(x@n)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@n[j])
                }
            }
        }, 
        "p" = {
            if (missing(j)) {
                return(x@p)
            } else {
                if (j>length(x@p)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@p[j])
                }
            }
        }, 
        "g" = {
            if (missing(j)) {
                return(x@g)
            } else {
                if (j>length(x@g)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@g[j])
                }
            }
        }, 
        "nItMC" = {
            if (missing(j)) {
                return(x@nItMC)
            } else {
                if (j>length(x@nItMC)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@nItMC[j])
                }
            }
        }, 
        "nItEM" = {
            if (missing(j)) {
                return(x@nItEM)
            } else {
                if (j>length(x@nItEM)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@nItEM[j])
                }
            }
        }, 
        "nBurn" = {
            if (missing(j)) {
                return(x@nBurn)
            } else {
                if (j>length(x@nBurn)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@nBurn[j])
                }
            }
        }, 
        "dp" = {
            if (missing(j)) {
                return(x@dp)
            } else {
                if (j>length(x@dp)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@dp[j])
                }
            }
        }, 
        "nsamp" = {
            if (missing(j)) {
                return(x@nsamp)
            } else {
                if (j>length(x@nsamp)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@nsamp[j])
                }
            }
        }, 
        "sparse" = {
            if (missing(j)) {
                return(x@sparse)
            } else {
                if (j>length(x@sparse)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@sparse[j])
                }
            }
        }, 
        "analysis" = {
            if (missing(j)) {
                return(x@analysis)
            } else {
                if (j>length(x@analysis)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@analysis[j])
                }
            }
        }, 
        "algorithm" = {
            if (missing(j)) {
                return(x@algorithm)
            } else {
                if (j>length(x@algorithm)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@algorithm[j])
                }
            }
        }, 
        "initialized" = {
            if (missing(j)) {
                return(x@initialized)
            } else {
                if (j>length(x@initialized)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@initialized[j])
                }
            }
        }, 
        "maxit" = {
            if (missing(j)) {
                return(x@maxit)
            } else {
                if (j>length(x@maxit)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@maxit[j])
                }
            }
        }, 
        "tol" = {
            if (missing(j)) {
                return(x@tol)
            } else {
                if (j>length(x@tol)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@tol[j])
                }
            }
        }, 
        "seed" = {
            if (missing(j)) {
                return(x@seed)
            } else {
                if (j>length(x@seed)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@seed[j])
                }
            }
        }, 
        "b" = {
            if (missing(j)) {
                return(x@b)
            } else {
                if (j>length(x@b)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
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
                    stop("[Clere:get] indice out of limits", call. = FALSE)
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
                    stop("[Clere:get] indice out of limits", call. = FALSE)
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
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@gamma2[j])
                }
            }
        }, 
        "intercept" = {
            if (missing(j)) {
                return(x@intercept)
            } else {
                if (j>length(x@intercept)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@intercept[j])
                }
            }
        }, 
        "likelihood" = {
            if (missing(j)) {
                return(x@likelihood)
            } else {
                if (j>length(x@likelihood)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
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
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@entropy[j])
                }
            }
        }, 
        "P" = {return(x@P)}, 
        "theta" = {return(x@theta)}, 
        "Zw" = {return(x@Zw)}, 
        "Bw" = {return(x@Bw)}, 
        "Z0" = {
            if (missing(j)) {
                return(x@Z0)
            } else {
                if (j>length(x@Z0)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@Z0[j])
                }
            }
        }, 
        "message" = {
            if (missing(j)) {
                return(x@message)
            } else {
                if (j>length(x@message)) {
                    stop("[Clere:get] indice out of limits", call. = FALSE)
                } else {
                    return(x@message[j])
                }
            }
        }, 
        stop("[Clere:get] ", i, " is not a \"Clere\" slot", call. = FALSE)
    )
})



### Clusters ###
setGeneric(name = "clusters", def = function(object, threshold = NULL, ...) {standardGeneric("clusters")})
setMethod(f = "clusters", signature = "Clere", definition = function(object, threshold = NULL, ...) {
    if(is.null(threshold)) {
      return(apply(object@P, 1, which.max))
    }else{
      return(sapply(1:object@p,function(j){
        counts = sum(object@P[j,]>=threshold)
        if( counts>0 ){
          if(counts>1){
            warning("[Clere:clusters] Variable ",j," could have been assigned to other groups using ",threshold," as threshold.\nYou may consider using a larger threshold or set option \" threshold = NULL\".",call.=TRUE)
          }
          return( which.max(object@P[j,])[1] )
        }else{
          return(NA)
        }
      }))
      ## return(apply(object@P, 1, function(x) ifelse(sum(x>=threshold)>0,which.max(x),NA)))
    }
})


### Predict ###
setMethod(f = "predict", signature = "Clere", definition = function(object, newx, ...) {
    if (class(newx) == "matrix") {
        if (ncol(newx) == object@p) {
          return(object@intercept + rowMeans(newx%*%object@Bw,na.rm=TRUE) )
        } else {
            stop(paste0("[Clere:predict] \"newx\" must be a matrix with ", object@p, " columns"), call. = FALSE)
        }
    } else {
        stop("[Clere:predict] \"newx\" is not a matrix", call. = FALSE)
    }
})


### Summary ###
setMethod(f = "summary", signature = "Clere", definition = function(object, ...) {
    if (missing(object)) {
        stop("[Clere:summary] \"object\" is missing", call. = FALSE)
    } else {}
    nbVarGroups = length(object@b)
    K = 2*(nbVarGroups+1)
    nd = 4
    sep = "\t"
    summaryClere <- new("sClere", 
        analysis = object@analysis, 
        g = object@g,
        nbVarGroups = length(object@b), 
        algorithm = object@algorithm, 
        intercept = object@intercept, 
        b = object@b, 
        pi = object@pi, 
        sigma2 = object@sigma2, 
        gamma2 = object@gamma2, 
        likelihood = object@likelihood, 
        entropy = object@entropy, 
        AIC = -2*object@likelihood + 2*K, 
        BIC = -2*object@likelihood + K*log(object@n), 
        ICL = -2*object@likelihood + K*log(object@n)+object@entropy
    )
    return(summaryClere)
})


### Plot ###
# setGeneric(name = "plot", def = function(x, ...) {standardGeneric("plot")})
setMethod(f = "plot", signature = "Clere", definition = function(x, ...) {
    if (nrow(x@theta) >= 2) {
        snbRow <- seq(nrow(x@theta))
        K <- length(x@b)
        sK <- seq(K)
        op <- par(mfrow = c(2, 2))
        xlabname = paste(x@algorithm, "iterations")
        
        ### First plot
        matplot(snbRow, x@theta[, 1+(sK)], type = "l", ylab = "The b's", lty = 1, xlab = xlabname)
        abline(v = x@nBurn, col = "grey")
        
        ### Second plot
        matplot(snbRow, x@theta[, 1+K+(sK)], type = "l", ylab = "The pi's", lty = 1, xlab = xlabname)
        abline(v = x@nBurn, col = "grey")
        
        ### Third plot
        matplot(snbRow, x@theta[, c(2*K+2, 2*K+3)], type = "l", ylab = "sigma^2 and gamma^2", lty = 1, xlab = xlabname, 
            ylim = c(0.5*min(x@theta[, c(2*K+2,2*K+3)], na.rm = TRUE), 1.5*max(x@theta[, c(2*K+2,2*K+3)], na.rm = TRUE)))
        legend("topright", legend = c("sigma2", "gamma2"), box.lty = 0, lty = 1, col = 1:2)
        abline(v = x@nBurn, col = "grey")
        
        ### Fourth plot
        plot(snbRow, x@theta[, 1], type = "l", ylab = "Intercept", lty = 1, xlab = xlabname)
        abline(v = x@nBurn, col = "grey")
        
        par(op)
    } else {
        stop("[Clere:plot] SEM iterations is below 2", call. = FALSE)
    }
})
