setMethod("sim", signature(object = "lm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    sigma.hat <- summ$sigma
    beta.hat <- coef[,1,drop = FALSE]
    V.beta <- summ$cov.unscaled
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    sigma <- rep (NA, n.sims)
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, rownames(beta.hat))
    for (s in 1:n.sims){
      sigma[s] <- sigma.hat*sqrt((n-k)/rchisq(1,n-k))
      beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta*sigma[s]^2)
    }

    ans <- new("sim",
                coef = beta,
                sigma = sigma)
    return (ans)
    }
)



setMethod("sim", signature(object = "glm"),
    function(object, n.sims=100)
    {
    object.class <- class(object)[[1]]
    summ <- summary (object, correlation=TRUE, dispersion = object$dispersion)
    coef <- summ$coef[,1:2,drop=FALSE]
    dimnames(coef)[[2]] <- c("coef.est","coef.sd")
    beta.hat <- coef[,1,drop=FALSE]
    sd.beta <- coef[,2,drop=FALSE]
    corr.beta <- summ$corr
    n <- summ$df[1] + summ$df[2]
    k <- summ$df[1]
    V.beta <- corr.beta * array(sd.beta,c(k,k)) * t(array(sd.beta,c(k,k)))
    beta <- array (NA, c(n.sims,k))
    dimnames(beta) <- list (NULL, dimnames(beta.hat)[[1]])
    for (s in 1:n.sims){
      beta[s,] <- MASS::mvrnorm (1, beta.hat, V.beta)
    }
    # Added by Masanao
    beta2 <- array (0, c(n.sims,length(coefficients(object))))
    dimnames(beta2) <- list (NULL, names(coefficients(object)))
    beta2[,dimnames(beta2)[[2]]%in%dimnames(beta)[[2]]] <- beta
    # Added by Masanao
    sigma <- rep (sqrt(summ$dispersion), n.sims)

    ans <- new("sim",
                coef = beta2,
                sigma = sigma)
    return(ans)
    }
)





setMethod("sim", signature(object = "polr"),
    function(object, n.sims=100){
  x <- as.matrix(model.matrix(object))
  coefs <- coef(object)
  k <- length(coefs)
  zeta <- object$zeta
  Sigma <- vcov(object)

  if(n.sims==1){
    parameters <- t(MASS::mvrnorm(n.sims, c(coefs, zeta), Sigma))
  }else{
    parameters <- MASS::mvrnorm(n.sims, c(coefs, zeta), Sigma)
  }
  ans <- new("sim.polr",
              coef = parameters[,1:k,drop=FALSE],
              zeta = parameters[,-(1:k),drop=FALSE])
  return(ans)
})



#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100)
#    {
#    #object <- summary(object)
##    if (lapply(object@bVar,sum)<=0|sum(unlist(lapply(object@bVar, is.na)))>0){
##        object@call$control <- list(usePQL=TRUE)
##        object <- lmer(object@call$formula)
#    #}
#    #sc <- attr (VarCorr (object), "sc")
#    # simulate unmodeled coefficients
#
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "unmodeled"
#    }
#    # simulate coefficients within groups
#    #coef <- ranef (object)
#    #estimate <- ranef(object, postVar=TRUE)
#    #vars <- object@bVar
#    #beta.bygroup <- vars
#
#    sc <- attr (VarCorr (object), "sc")
#    coef <- ranef(object, postVar=TRUE)
#    beta.bygroup <- c(coef)
#    n.groupings <- length (coef)
#    for (m in 1:n.groupings){
#      #vars.m <- vars[[m]]
#      vars.m <- attr (coef[[m]], "postVar")
#      K <- dim(vars.m)[1]
#      J <- dim(vars.m)[3]
#      beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#      bhat <- coef[[m]]
#      for (j in 1:J){
#        V.beta <- untriangle(vars.m[,,j])#*sc^2
#        beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#      }
#      dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#    }
#    betas <- c (beta.unmodeled, beta.bygroup)
#    return (betas)
#    }
#)

#setMethod("sim", signature(object = "mer"),
#    function(object, n.sims=100, ranef=TRUE)
#    {
#    # simulate unmodeled coefficients
#    fcoef <- fixef(object)
#    corF <- vcov(object)@factors$correlation
#    se.unmodeled <- corF@sd
#    V.beta <- (se.unmodeled %o% se.unmodeled) * as.matrix(corF)
#    beta.unmodeled <- NULL
#    if (length (fcoef) > 0){
#      beta.unmodeled[[1]] <- mvrnorm (n.sims, fcoef, V.beta)
#      names (beta.unmodeled) <- "fixef"#"unmodeled"
#      coef <- beta.unmodeled
#    }
#    if(ranef){
#      # simulate coefficients within groups
#      sc <- attr (VarCorr (object), "sc")  # scale
#      #coef <- ranef (object)
#      #estimate <- ranef(object, postVar=TRUE)
#      coef <- ranef(object, postVar=TRUE)
#      beta.bygroup <- coef
#      n.groupings <- length (coef)
#      for (m in 1:n.groupings){
#        bhat <- as.matrix(coef[[m]]) # to suit the use of mvrnorm
#        vars.m <- attr (coef[[m]], "postVar")
#        K <- dim(vars.m)[1]
#        J <- dim(vars.m)[3]
#        beta.bygroup[[m]] <- array (NA, c(n.sims, J, K))
#        for (j in 1:J){
#          V.beta <- .untriangle(vars.m[,,j])#*sc^2
#          beta.bygroup[[m]][,j,] <- mvrnorm (n.sims, bhat[j,], V.beta)
#        }
#        dimnames (beta.bygroup[[m]]) <- c (list(NULL), dimnames(bhat))
#      }
#      coef <- c (beta.unmodeled, beta.bygroup)
#      }
#    return (coef)
#    }
#)
