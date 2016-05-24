### This file contains functions to call phyclust_em_step, phyclust_e_step, and
### phyclust_m_step in C.

### EM step.
phyclust.em.step <- function(X, ret.phyclust = NULL,
    K = NULL, Eta = NULL, Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL){
  if(is.null(ret.phyclust)){
    if(is.null(K) || is.null(Eta) || is.null(Mu) || is.null(Tt) ||
       is.null(substitution.model) ||
       is.null(identifier) || is.null(code.type)){
      stop("The parameters are not specified correctly.")
    } else{
      ret.phyclust <- list(K = K, Eta = Eta, Mu = Mu,
                           QA = list(pi = pi, kappa = kappa, Tt = Tt,
                                     identifier = identifier),
                           substitution.model = substitution.model,
                           code.type = code.type)
    }
  } else{
    if(class(ret.phyclust) != "phyclust"){
      stop("The ret.phyclust should be in a phyclust class.")
    }
  }

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, nrow(X), ret.phyclust$K, TRUE)

  ret <- .Call("R_phyclust_em_step",
               as.integer(nrow(X)),
               as.integer(ncol(X)),
               as.integer(t(X)),
               as.integer(ret.phyclust$K),
               as.double(ret.phyclust$Eta),
               as.integer(t(ret.phyclust$Mu)),
               as.double(vect),
               as.integer(which(ret.phyclust$substitution.model ==
                                as.character(.substitution.model$model)) - 1),
               as.integer(which(ret.phyclust$QA$identifier == .identifier) - 1),
               as.integer(which(ret.phyclust$code.type == .code.type) - 1),
               label,
               PACKAGE = "phyclust")

  ret$Z.normalized <- ret$bic <- ret$aic <- ret$icl <-
    ret$class.id <- ret$n.class <- NULL
  ret$substitution.model <- ret.phyclust$substitution.model
  ret$QA$identifier <- ret.phyclust$QA$identifier
  ret$code.type <- ret.phyclust$code.type

  ret <- translate.ret(ret)
  class(ret) <- "phyclust"
  ret
} # End of phyclust.em.step().


### E-step: return a matrix Z.normalized with dim=NxK.
### Z.state = 1 return Z.normalized
###           0 return logPt
phyclust.e.step <- function(X, ret.phyclust = NULL,
    K = NULL, Eta = NULL, Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    Z.state = TRUE, label = NULL){
  if(is.null(ret.phyclust)){
    if(is.null(K) || is.null(Eta) || is.null(Mu) || is.null(Tt) ||
       is.null(substitution.model) ||
       is.null(identifier) || is.null(code.type)){
      stop("The parameters are not specified correctly.")
    } else{
      ret.phyclust <- list(K = K, Eta = Eta, Mu = Mu,
                           QA = list(pi = pi, kappa = kappa, Tt = Tt,
                                     identifier = identifier),
                           substitution.model = substitution.model,
                           code.type = code.type)
    }
  } else{
    if(class(ret.phyclust) != "phyclust"){
      stop("The ret.phyclust should be in a phyclust class.")
    }
  }

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, nrow(X), ret.phyclust$K, TRUE)

  ret <- .Call("R_phyclust_e_step",
               as.integer(nrow(X)),
               as.integer(ncol(X)),
               as.integer(t(X)),
               as.integer(ret.phyclust$K),
               as.double(ret.phyclust$Eta),
               as.integer(t(ret.phyclust$Mu)),
               as.double(vect),
               as.integer(which(ret.phyclust$substitution.model ==
                                as.character(.substitution.model$model)) - 1),
               as.integer(which(ret.phyclust$QA$identifier == .identifier) - 1),
               as.integer(which(ret.phyclust$code.type == .code.type) - 1),
               as.integer(Z.state),
               label,
               PACKAGE = "phyclust")
  ret <- matrix(ret, nrow = nrow(X), byrow = TRUE)

  ret
} # End of phyclust.e.step().


### M-step: return a object with phyclust class.
phyclust.m.step <- function(X, ret.phyclust = NULL,
    K = NULL, pi = NULL, kappa = NULL, Tt = NULL, Z.normalized = NULL,
    substitution.model = NULL, identifier = NULL, code.type = NULL,
    label = NULL){
  if(is.null(ret.phyclust)){
    if(is.null(K) || is.null(Tt) ||
       is.null(Z.normalized) || is.null(substitution.model) ||
       is.null(identifier) || is.null(code.type)){
      stop("The parameters are not specified correctly.")
    } else{
      ret.phyclust <- list(K = K,
                           QA = list(pi = pi, kappa = kappa, Tt = Tt,
                                     identifier = identifier),
                           Z.normalized = Z.normalized,
                           substitution.model = substitution.model,
                           code.type = code.type)
    }
  } else{
    if(class(ret.phyclust) != "phyclust"){
      stop("The ret.phyclust should be in a phyclust class.")
    }
  }

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, nrow(X), ret.phyclust$K, TRUE)

  ret <- .Call("R_phyclust_m_step",
               as.integer(nrow(X)),
               as.integer(ncol(X)),
               as.integer(t(X)),
               as.integer(ret.phyclust$K),
               as.double(vect),
               as.double(t(ret.phyclust$Z.normalized)),
               as.integer(which(ret.phyclust$substitution.model ==
                                as.character(.substitution.model$model)) - 1),
               as.integer(which(ret.phyclust$QA$identifier == .identifier) - 1),
               as.integer(which(ret.phyclust$code.type == .code.type) - 1),
               label,
               PACKAGE = "phyclust")

  ret$Z.normalized <- ret$bic <- ret$aic <- ret$icl <-
    ret$class.id <- ret$n.class <- NULL
  ret$substitution.model <- ret.phyclust$substitution.model
  ret$QA$identifier <- ret.phyclust$QA$identifier
  ret$code.type <- ret.phyclust$code.type

  ret <- translate.ret(ret)
  class(ret) <- "phyclust"
  ret
} # End of phyclust.m.step().

