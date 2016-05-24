setMethod("se.coef", signature(object = "lm"),
    function(object)
    {
    object.class <- class(object)[[1]]
    sqrt (diag(vcov(object)))
    }
)


setMethod("se.coef", signature(object = "glm"),
    function(object)
    {
    object.class <- class(object)[[1]]
    sqrt (diag(vcov(object)))
    }
)

#setMethod("se.coef", signature(object = "mer"),
#    function(object)
#    {
#    #    if (sum(unlist(lapply(object@bVar, is.na)))>0){
##        object@call$control <- list(usePQL=TRUE)
##        object <- lmer(object@call$formula)
##    }
#    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
#    fcoef <- fixef(object)
#    #sc <- attr (VarCorr (object), "sc")
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- NULL
#    se.unmodeled[[1]] <- corF@sd
#    names (se.unmodeled) <- "unmodeled"
#
#    #coef <- ranef (object)
#    #estimate <- ranef(object, postVar=TRUE)
#    coef <- ranef(object, postVar=TRUE)
#    se.bygroup <- coef #ranef( object, postVar = TRUE )
#    n.groupings <- length (coef)
#    
#    for (m in 1:n.groupings){
#      vars.m <- attr (coef[[m]], "postVar")
#      K <- dim(vars.m)[1]
#      J <- dim(vars.m)[3]
#      se.bygroup[[m]] <- array (NA, c(J,K))
#      for (j in 1:J){
#        se.bygroup[[m]][j,] <- sqrt(diag(as.matrix(vars.m[,,j])))
#      }
##      se.bygroup[[m]] <- se.bygroup[[m]]*sc
#      names.full <- dimnames (ranef(object)[[m]])
#      dimnames (se.bygroup[[m]]) <- list (names.full[[1]],
#                            names.full[[2]])
#    }
#    #names(se.bygroup) <- names(ngrps)
#    ses <- c (se.unmodeled, se.bygroup)
#    return (ses)
#    }
#)

setMethod("se.coef", signature(object = "merMod"),
    function(object)
    {
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    fcoef <- fixef(object)
    #sc <- attr (VarCorr (object), "sc")
    corF <- vcov(object)@factors$correlation
    se.unmodeled <- NULL
    se.unmodeled[[1]] <- corF@sd
    names (se.unmodeled) <- "fixef"#"unmodeled"

    #coef <- ranef (object)
    #estimate <- ranef(object, postVar=TRUE)
    coef <- ranef(object, condVar=TRUE)
    se.bygroup <- coef #ranef( object, postVar = TRUE )
    n.groupings <- length (coef)
    
    for (m in 1:n.groupings){
      vars.m <- attr (coef[[m]], "postVar")
      K <- dim(vars.m)[1]
      J <- dim(vars.m)[3]
      se.bygroup[[m]] <- array (NA, c(J,K))
      for (j in 1:J){
        se.bygroup[[m]][j,] <- sqrt(diag(as.matrix(vars.m[,,j])))
      }
#      se.bygroup[[m]] <- se.bygroup[[m]]*sc
      names.full <- dimnames (coef[[m]])
      dimnames (se.bygroup[[m]]) <- list (names.full[[1]],
                            names.full[[2]])
    }
    #names(se.bygroup) <- names(ngrps)
    ses <- c (se.unmodeled, se.bygroup)
    return (ses)
    }
)



se.fixef <- function (object){
  #object <- summary (object)
  fcoef.name <- names(fixef(object))
  corF <- vcov(object)@factors$correlation
  ses <- corF@sd
  names(ses) <- fcoef.name
  return (ses)
}

se.ranef <- function (object){
    #ngrps <- lapply(object@flist, function(x) length(levels(x)))
    se.bygroup <- ranef( object, condVar = TRUE )
    n.groupings<- length( se.bygroup )
    for( m in 1:n.groupings ) {
        vars.m <- attr( se.bygroup[[m]], "postVar" )
        K <- dim(vars.m)[1]
        J <- dim(vars.m)[3]
        names.full <- dimnames(se.bygroup[[m]])
        se.bygroup[[m]] <- array(NA, c(J, K))
        for (j in 1:J) {
            se.bygroup[[m]][j, ] <- sqrt(diag(as.matrix(vars.m[, , j])))
        }       
        dimnames(se.bygroup[[m]]) <- list(names.full[[1]], names.full[[2]])
    }
    return(se.bygroup)
}
