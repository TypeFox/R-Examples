predict.DirichletRegModel <- function(object, newdata, mu = TRUE, alpha = FALSE, phi = FALSE, ...){

  if(missing(newdata)) return(fitted(object, mu, alpha, phi))

  repar <- object$parametrization == "alternative"
  dims  <- ncol(object$Y)

  model_formula <- object$mf_formula
  model_formula$formula <- as.Formula(deparse(model_formula$formula)) # circumvent problems with Formula objects
  model_formula$data <- as.name("newdata")
  model_formula$lhs <- 0

  # include intercept for repar-model w/o phi block
  if(repar && (length(model_formula$formula)[2L] == 1L)){
    model_formula$formula <- as.Formula(paste0(deparse(model_formula$formula), " | 1"))
  }

  # fill all rhs blocks if only one is specified
  if(!repar && (length(model_formula$formula)[2L] == 1L)){
    model_formula$formula <- as.Formula( paste0(deparse(model_formula$formula),
      " | ", paste0(rep(deparse(model_formula$formula[[3]]), dims - 1L), collapse=" | ")) )
  }

  model_formula[["drop.unused.levels"]] <- FALSE
  mf <- eval(model_formula)

  if(!repar){ ## COMMON
    X <- lapply(seq_len(dims), function(i){ model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=i)), mf) })
    Z <- NULL
  } else { ## ALTERNATIVE
    X <- model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=1L)), mf)
    Z <- model.matrix(Formula(terms(model_formula$formula, data=newdata, rhs=2L)), mf)
  }

  cc <- coef(object)

  if(repar){
    base <- object$base
    cc[[1L]] <- split(unlist(cc[[1L]]), factor(seq_len(dims))[rep(seq_len(dims)[-base], each=ncol(X))])
    cc[[2L]] <- unlist(cc[[2L]])

    ETA <- matrix(0, nrow = nrow(newdata), ncol = dims)
    for(i in seq_len(dims)[-base]){
      ETA[,i] <- X %*% cc[[1]][[i]]
    }
    MU    <- exp(ETA)/rowSums(exp(ETA))
    PHI   <- exp(Z %*% cc[[2L]])
    ALPHA <- MU * as.numeric(PHI)
  } else {
    ALPHA <- matrix(0, nrow = nrow(newdata), ncol = dims)
    for(i in seq_len(dims)){
      ALPHA[,i] <- exp(X[[i]] %*% cc[[i]])
    }
    PHI <- rowSums(ALPHA)
    MU <- ALPHA/PHI
  }

  if(!any(mu || alpha || phi)) stop("Either mu, alpha or phi has to be requested.")

  if(sum(mu + alpha + phi) == 1){
    if(mu)    return(MU)
    if(alpha) return(ALPHA)
    if(phi)   return(PHI)
  } else {
    res <- list()
    if(mu)    res[["mu"]]    <- MU
    if(alpha) res[["alpha"]] <- ALPHA
    if(phi)   res[["phi"]]   <- PHI

    return(res)
  }

}
