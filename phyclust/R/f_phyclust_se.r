# This file contains functions to call phyclust in C.

phyclust.se <- function(X, K, EMC = .EMC, manual.id = NULL, byrow = TRUE){
  label <- NULL

  if(K <= 0){
    stop("K > 0")
  }

  if(! is.null(manual.id)){
    if(! is.null(label)){
      stop("The manual.id is only for unsupervised clustering.")
    }
    if(any(manual.id <= 0 | manual.id > K)){
      stop("The manual.id is not correct.")
    }
    manual.id <- manual.id - 1
  } else{
    if(EMC$init.method == "manualMu"){
      stop("The manual.id is missing.")
    }
  }

  if(byrow){
    X <- t(X)
  }
  N.X.org <- ncol(X)
  L <- nrow(X)

  EMC <- check.EMC(EMC)
  EMC <- translate.EMC.se(EMC)
  label <- check.label(label, N.X.org, K, byrow)

  ret <- .Call("R_phyclust_se",
               as.integer(N.X.org),
               as.integer(L),
               as.integer(K),
               as.integer(X),
               EMC,
               as.integer(manual.id),
               label,
               PACKAGE = "phyclust")

  if(!is.finite(ret$logL)){
    stop("The logL is not finite.\n")
  }

  ret$class.id <- ret$class.id + 1
  ret <- translate.ret.se(ret, EMC)
  class(ret) <- "phyclust"
  ret
} # End of phyclust().


### For internal used.
translate.ret.se <- function(ret, EMC = NULL){
  ret$Mu <- matrix(ret$Mu, nrow = ret$K, byrow = TRUE)
  if(! is.null(ret$Z.normalized)){
    ret$Z.normalized <- matrix(ret$Z.normalized, nrow = ret$N.X.org, byrow = TRUE)
  }
  ret$QA$pi <- matrix(ret$QA$pi, nrow = ret$K, byrow = TRUE)
  ret$QA$kappa <- matrix(ret$QA$kappa, nrow = ret$K, byrow = TRUE)
  ret$QA$Tt <- matrix(ret$QA$Tt, nrow = ret$K, byrow = TRUE)

  if(!is.null(EMC)){
    ret$init.procedure <- .init.procedure[EMC$init.procedure + 1]
    ret$init.method <- .init.method[EMC$init.method + 1]
    ret$substitution.model <-
      as.character(.substitution.model$model[EMC$substitution.model + 1])
    ret$edist.model <- .edist.model[EMC$edist.model + 1]
    ret$QA$identifier <- .identifier[EMC$identifier + 1]
    ret$code.type <- .code.type[EMC$code.type + 1]
    ret$em.method <- .em.method[EMC$em.method + 1]
    ret$boundary.method <- .boundary.method[EMC$boundary.method + 1]
  }

  if(ret$QA$identifier == "EE"){
    ret$QA$pi <- matrix(ret$QA$pi[1,], nrow = 1)
    ret$QA$kappa <- ret$QA$kappa[1]
    ret$QA$Tt <- ret$QA$Tt[1]
  } else if(ret$QA$identifier == "EV"){
    ret$QA$pi <- matrix(ret$QA$pi[1,], nrow = 1)
    ret$QA$kappa <- ret$QA$kappa[1]
  } else if(ret$QA$identifier == "VE"){
    ret$QA$Tt <- ret$QA$Tt[1]
  }

  if(ret$code.type == "NUCLEOTIDE"){
    colnames(ret$QA$pi) <- as.character(.nucleotide$code[1:ncol(ret$QA$pi)])
  } else if(ret$code.type == "SNP"){
    colnames(ret$QA$pi) <- as.character(.snp$code[1:ncol(ret$QA$pi)])
  }
  rownames(ret$QA$pi) <- paste("k=", 1:nrow(ret$QA$pi), sep = "")

  if(ret$substitution.model %in% c("JC69", "K80", "SNP_JC69")){
    ret$QA$pi <- NULL
  }
  if(ret$substitution.model %in%
     c("JC69", "F81", "SNP_JC69", "SNP_F81", "E_F81", "E_SNP_F81")){
    ret$QA$kappa <- NULL
  }

  # ret$label.method <- .label.method[ret$label.method + 1]

  ### For se model.
  ret$SE$type <- as.logical(ret$SE$type)
  ret$SE$model <- .se.model[ret$SE$model + 1]
  ret$SE$f.err <- matrix(ret$SE$f.err, nrow = 4, byrow = TRUE)
  colnames(ret$SE$f.err) <- as.character(.nucleotide$code[1:ncol(ret$SE$f.err)])
  rownames(ret$SE$f.err) <- as.character(.nucleotide$code[1:nrow(ret$SE$f.err)])

  ret
} # End of translate.ret.se().

translate.EMC.se <- function(EMC){
  EMC$exhaust.iter <- as.integer(EMC$exhaust.iter)
  EMC$fixed.iter <- as.integer(EMC$fixed.iter)
  EMC$short.iter <- as.integer(EMC$short.iter)
  EMC$EM.iter <- as.integer(EMC$EM.iter)
  EMC$short.eps <- as.double(EMC$short.eps)
  EMC$EM.eps <- as.double(EMC$EM.eps)

  EMC$cm.reltol <- as.double(EMC$cm.reltol)
  EMC$cm.maxit <- as.integer(EMC$cm.maxit)

  EMC$nm.abstol.Mu.given.QA <- as.double(EMC$nm.abstol.Mu.given.QA)
  EMC$nm.abstol.QA.given.Mu <- as.double(EMC$nm.abstol.QA.given.Mu)
  EMC$nm.reltol.Mu.given.QA <- as.double(EMC$nm.reltol.Mu.given.QA)
  EMC$nm.reltol.QA.given.Mu <- as.double(EMC$nm.reltol.QA.given.Mu)
  EMC$nm.maxit.Mu.given.QA <- as.integer(EMC$nm.maxit.Mu.given.QA)
  EMC$nm.maxit.QA.given.Mu <- as.integer(EMC$nm.maxit.QA.given.Mu)
  EMC$est.non.seg.site <- as.integer(EMC$est.non.seg.site)

  EMC$min.n.class <- as.integer(EMC$min.n.class)

  if(EMC$init.procedure[1] %in% .init.procedure){
    EMC$init.procedure <- as.integer(which(EMC$init.procedure ==
                                           .init.procedure) - 1)
  } else{
    stop("The initial procedure is not found.")
  }

  if(EMC$init.method[1] %in% .init.method){
    EMC$init.method <- as.integer(which(EMC$init.method == .init.method) - 1)
  } else{
    stop("The initial method is not found.")
  }

  if(EMC$substitution.model[1] %in% as.character(.substitution.model$model)){
    EMC$substitution.model <-
      as.integer(which(EMC$substitution.model ==
                       as.character(.substitution.model$model)) - 1)
  } else{
    stop("The substitution model is not found.")
  }

  if(EMC$edist.model[1] %in% .edist.model){
    EMC$edist.model <- as.integer(which(EMC$edist.model == .edist.model) - 1)
  } else{
    stop("The edist model is not found.")
  }

  if(EMC$identifier[1] %in% .identifier){
    EMC$identifier <- as.integer(which(EMC$identifier == .identifier) - 1)
  } else{
    stop("The identifier is not found.")
  }

  if(EMC$code.type[1] %in% .code.type){
    EMC$code.type <- as.integer(which(EMC$code.type == .code.type) - 1)
  } else{
    stop("The code type is not found.")
  }

  if(EMC$em.method[1] %in% .em.method){
    EMC$em.method <- as.integer(which(EMC$em.method == .em.method) - 1)
  } else{
    stop("The em method is not found.")
  }

  if(EMC$boundary.method[1] %in% .boundary.method){
    EMC$boundary.method <- as.integer(which(EMC$boundary.method ==
                                            .boundary.method) - 1)
  } else{
    stop("The boundary method is not found.")
  }

  EMC$max.init.iter <- as.integer(EMC$max.init.iter)

  ### For se model.
  if(EMC$code.type == as.integer(which(.code.type == "NUCLEOTIDE") - 1)){
    EMC$se.type <- as.integer(1)
    EMC$se.model <- as.integer(which(EMC$se.model == .se.model) - 1)
    EMC$se.constant <- as.double(EMC$se.constant)
  } else{
    stop("The sequencing error model only supports NUCLEOTIDE data.")
  }

  EMC
} # End of translate.EMC.se().

