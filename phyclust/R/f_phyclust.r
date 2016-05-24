# This file contains functions to call phyclust in C.

phyclust <- function(X, K, EMC = .EMC, manual.id = NULL,
    label = NULL, byrow = TRUE){
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
  EMC <- translate.EMC(EMC)
  label <- check.label(label, N.X.org, K, byrow)

  ret <- .Call("R_phyclust",
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
  ret <- translate.ret(ret, EMC)
  class(ret) <- "phyclust"
  ret
} # End of phyclust().


### For internal used.
check.EMC <- function(EMC){
  if(EMC$init.method %in% c("NJ", "PAM", "manualMu") &&
     EMC$init.procedure != "exhaustEM"){
    my.cat("init procedure: ", EMC$init.procedure, " -> exhaustEM",
           " (", EMC$init.method, ")\n")
    EMC$init.procedure <- "exhaustEM"
    EMC$exhaust.iter <- 1
  }

  if((EMC$code.type == "SNP") &&
     (!EMC$substitution.model %in% c("SNP_JC69", "SNP_F81", "E_SNP_F81"))){
    my.cat("substitution model: ", EMC$substitution.model, " -> SNP_JC69",
           " (", EMC$code.type, ")\n")
    EMC$substitution.model <- "SNP_JC69"
  }

  if(EMC$substitution.model %in% c("SNP_JC69", "SNP_F81", "E_SNP_F81") &&
     ! (EMC$edist.model %in% c("D_HAMMING", "D_HAMMING_WOGAP"))){
    my.cat("edist model: ", EMC$edist.model, " -> D_HAMMING",
           " (", EMC$substitution.model, ")\n")
    EMC$edist.model <- "D_HAMMING"
  }

  EMC
} # End of check.EMC().

translate.ret <- function(ret, EMC = NULL){
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

  ret$label.method <- .label.method[ret$label.method + 1]

  ret
} # End of translate.ret().

translate.EMC <- function(EMC){
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

  EMC
} # End of translate.EMC().


### For internal used in semi-supervised clustering.
check.label <- function(label, N.X.org, K, byrow){
  ### label$label.method = 0 for NONE, 1 for SIMPLE, and 2 for COMPLEX.
  if(is.null(label)){
    ### Do nothing for unsupervised.
    label <- list(label.method = as.integer(which(.label.method == "NONE") - 1),
               semi = NULL, index = NULL, prob = NULL)
  } else{
    ### Check for semi-supervised.
    if(is.vector(label, mode = "numeric")){
      ### A vector of label is for semi-supervised clustering.
      label <- list(label.method = as.integer(which(.label.method == "SEMI") - 1),
                 semi = label, index = NULL, prob = NULL)

      if((length(label$semi) != N.X.org) ||
         (length(unique(label$semi)) != (max(label$semi) + 1)) ||
         any(label$semi < 0 | label$semi > K)){
        stop("The label$semi is not correct.")
      }

      label$index <- which(label$semi != 0)
      label$prob <- matrix(0, nrow = K, ncol = length(label$index))
      for(i in 1:length(label$index)){
        label$prob[label$semi[label$index[i]], i] <- 1
      }
    } else{
      ### A data.frame of label is for semi-supervised clustering.
      if(! is.vector(label, mode = "list")){
        stop("The label should be a list.")
      }

      label$label.method <- as.integer(which(.label.method == "GENERAL") - 1)

      if(is.null(label$index) || is.null(label$prob)){
        stop("The label$index and label$prob are required.")
      }
      if(length(label$index) > N.X.org ||
         any(label$index < 1 | label$index > N.X.org)){
        stop("The label$index is not correct.")
      }

      if(byrow){
        label$prob <- t(label$prob)
      }
      label$semi <- apply(label$prob, 2, which.max)

      if((nrow(label$prob) != K) ||
         (ncol(label$prob) != length(label$index)) ||
         any(label$prob < 0 | label$prob > 1) ||
         any(colSums(label$prob) != 1)){
        stop("The label$prob is not correct.")
      }
    }

    label$semi <- as.integer(label$semi[label$semi > 0] - 1)
    label$index <- as.integer(label$index - 1)
  }

  label
} # End of check.label().

