# This file contains functions to call phyclust_logL in C.

phyclust.logL <- function(X, ret.phyclust = NULL,
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

  if(ret.phyclust$substitution.model == "E_F81"){
    ret.phyclust$substitution.model <- "F81"
  } else if(ret.phyclust$substitution.model == "E_HKY85"){
    ret.phyclust$substitution.model <- "HKY85"
  } else if(ret.phyclust$substitution.model == "E_SNP_F81"){
    ret.phyclust$substitution.model <- "SNP_F81"
  }

  vect <- convert.QA.to.vect(ret.phyclust)
  label <- check.label(label, nrow(X), ret.phyclust$K, TRUE)

  ret <- .Call("R_phyclust_logL",
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

  ret
} # End of phyclust.logL().


### For internal used.
# Storage of vect for EE:
# 	0 to (n_param - 1):				for k = 0, 1, 2,..., K-1.
# Storage of vect for EV:
# 	0 to (n_param - 2):				for k = 0, 1, 2,..., K-1.
#	(n_param - 1) to (n_param - 1 + K):		for Tt[0] to Tt[K-1].
# Storage of vect for VE:
# 	0 to (n_param - 2):				for k = 0,
#	(n_param - 1) to (2*(n_param - 1) - 1):		for k = 1,
#	(2*(n_param - 1)) to (3*(n_param - 1) - 1):	for k = 2,
#	... until k = K-1.
#	((K-1)*(n_param - 1))				for Tt.
# Storage of vect for VV:
#	0 to (n_param - 1):				for k = 0,
#	(n_param) to (2*n_param - 1):			for k = 1,
#	(2*n_param) to (3*n_param - 1):			for k = 2,
# 	... until k = K-1.
convert.QA.to.vect <- function(ret.phyclust){
  if(ret.phyclust$code.type == "NUCLEOTIDE"){
    ncode <- 3
  } else if(ret.phyclust$code.type == "SNP"){
    ncode <- 1
  }
  model.with.pi <- c("F81", "HKY85", "SNP_F81")
  model.with.kappa <- c("K80", "F81", "HKY85", "SNP_F81")

  K <- ret.phyclust$K
  vect <- NULL

  if(ret.phyclust$QA$identifier == "EE"){
    if(ret.phyclust$substitution.model %in% model.with.pi){
      vect <- c(vect, ret.phyclust$QA$pi[1:ncode])
    }
    if(ret.phyclust$substitution.model %in% model.with.kappa){
      vect <- c(vect, ret.phyclust$QA$kappa)
    }
    vect <- c(vect, ret.phyclust$QA$Tt)
  } else if(ret.phyclust$QA$identifier == "EV"){
    if(ret.phyclust$substitution.model %in% model.with.pi){
      vect <- c(vect, ret.phyclust$QA$pi[1:ncode])
    }
    if(ret.phyclust$substitution.model %in% model.with.kappa){
      vect <- c(vect, ret.phyclust$QA$kappa)
    }
    vect <- c(vect, ret.phyclust$QA$Tt)
  } else if(ret.phyclust$QA$identifier == "VE"){
    for(i in 1:K){
      if(ret.phyclust$substitution.model %in% model.with.pi){
        vect <- c(vect, ret.phyclust$QA$pi[i, 1:ncode])
      }
      if(ret.phyclust$substitution.model %in% model.with.kappa){
        vect <- c(vect, ret.phyclust$QA$kappa[i])
      }
    }
    vect <- c(vect, ret.phyclust$QA$Tt)
  } else{
    for(i in 1:K){
      if(ret.phyclust$substitution.model %in% model.with.pi){
        vect <- c(vect, ret.phyclust$QA$pi[i, 1:ncode])
      }
      if(ret.phyclust$substitution.model %in% model.with.kappa){
        vect <- c(vect, ret.phyclust$QA$kappa[i])
      }
      vect <- c(vect, ret.phyclust$QA$Tt[i])
    }
  }

  vect
} # End of convert.QA.to.vect().

