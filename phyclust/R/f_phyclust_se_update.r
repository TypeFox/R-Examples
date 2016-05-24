### This file contains functions to call phyclust_se_update in C.

### EM step.
phyclust.se.update <- function(X, EMC = .EMC, ret.phyclust = NULL,
    K = NULL, Eta = NULL, Mu = NULL, pi = NULL, kappa = NULL, Tt = NULL,
    byrow = TRUE){

  label <- NULL

  if(is.null(ret.phyclust)){
    if(is.null(K) || is.null(Eta) || is.null(Mu) || is.null(Tt)){
      stop("The parameters are not specified correctly.")
    } else{
      ret.phyclust <- list(K = K, Eta = Eta, Mu = Mu,
                           QA = list(pi = pi, kappa = kappa, Tt = Tt,
                                     identifier = EMC$identifier),
                           substitution.model = EMC$substitution.model,
                           code.type = EMC$code.type)
    }
  } else{
    if(class(ret.phyclust) != "phyclust"){
      stop("The ret.phyclust should be in a phyclust class.")
    }
  }

  if(ret.phyclust$code.type != "NUCLEOTIDE"){
    stop("The sequencing error model only supports NUCLEOTIDE data.")
  }

  K <- ret.phyclust$K
  if(byrow){
    X <- t(X)
  }
  N.X.org <- ncol(X)
  L <- nrow(X)

  EMC$init.procedure <- ret.phyclust$init.procedure
  EMC$init.method <- ret.phyclust$init.method
  EMC$substitution.model <- ret.phyclust$substitution.model
  EMC$edist.model <- ret.phyclust$edist.model
  EMC$identifier <- ret.phyclust$QA$identifier
  EMC$code.type <- ret.phyclust$code.type
  EMC$em.method <- "EM"
  EMC$boundary.method <- ret.phyclust$boundary.method
  EMC$se.type <- TRUE
  EMC$se.model <- "CONVOLUTION"

  EMC <- check.EMC(EMC)
  EMC <- translate.EMC.se(EMC)

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, N.X.org, K, byrow)

  ret <- .Call("R_phyclust_se_update",
               as.integer(N.X.org),
               as.integer(L),
               as.integer(X),
               EMC,
               as.integer(K),
               as.double(ret.phyclust$Eta),
               as.integer(t(ret.phyclust$Mu)),
               as.double(vect),
               label,
               PACKAGE = "phyclust")

  if(!is.finite(ret$logL)){
    stop("The logL is not finite.\n")
  }

  ret$class.id <- ret$class.id + 1
  ret <- translate.ret.se(ret, EMC)
  class(ret) <- "phyclust"
  ret
} # End of phyclust.se.update().

