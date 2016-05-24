print.markophylo <- function(x, ...) {
  cat("Call:\n")
  print.default(x$call)
  cat("\n", x$taxa, "taxa with", x$phyl, "discrete character patterns were fitted.\n")
  cat("-----------------------------------\n")
  cat("Groups of nodes with the same rates:\n")
  print.default(x$bg)
  cat("-----------------------------------\n")
  cat("Parameter Index Matrix\n")
  mm_al <- x$modelmat
  rownames(mm_al) <- colnames(mm_al) <- x$alphabet
  print(mm_al)
  cat("-----------------------------------\n")
  for (i in x$modelnames) {
    cat(i, "\n")
    if (x$results[[i]]$convergence == 0) {
      print.default(x$results[[i]]$parsep)
      cat("-----------------------------------\n")
      if (x$reversible) {
        cat("Reversible option used. The root probabilities are estimated using maximum likelihood (see above for estimates and standard errors), and the rate matrix supplied to the function is symmetric (see above for estimates and standard errors). The estimated rate matrix can be written as the product of the symmetric matrix multiplied by a diagonal matrix (consisting of the estimated root probabilities). These rate estimates are:\n")
        mm_al2 <- x$results[[i]]$modelmat2
        rownames(mm_al2) <- colnames(mm_al2) <- x$alphabet
        print(mm_al2)
        print.default(x$results[[i]]$revcombined)
      }
      cat("Log-likelihood for model", i, ":", -x$results[[i]]$objective, "\n")
      cat("AIC           for model", i, ":", x$results[[i]]$AIC, "\n")
      cat("BIC           for model", i, ":", x$results[[i]]$BIC, "\n")
      cat("-----------------------------------\n")
    } else {
      cat("Convergence not achieved.\nIncrease number of iterations and/or function calls. See control arguments for the optimization method used.\n")
      cat("-----------------------------------\n")
    }
  }
  if (!all(x$conv == 0)) 
    cat("Not all models converged. \n")
  cat("Time taken:", x$time[3], "seconds.\n")
}
# calibrate_ratemat <- function(x, Q, iroot){
#   #Outputs the average number of transition events per site per evolutionary time unit.
#   #x takes ...$results$wop$parsep$rates[,,1]
#   #OR for reversible x takes x$results$wop$revcombined$rates[,,1]
#   #Q takes ...$modelmat or ...$results$wop$modelmat2
#   #iroot takes ...$results$wop$parsep$rates[,,1]
#   #Should probably output iroot no matter if it is equal, stationary, or rpvec.
#   #Already output if maxlik
#   Q_out <- Q
#   for (gf in 1:length(x)) {
#     Q_out[Q == gf] <- x[gf]
#   }
#   diag(Q_out) <- -rowSums(Q_out, na.rm = TRUE)
#   return(-sum(diag(Q_out)*iroot))
# }
plottree <- function(x, colors = NULL, ...) {
  if (is.null(colors)) 
    colors <- grDevices::rainbow(length(x$bg))
  colo <- NULL
  for (i in 1:nrow(x$tree$edge)) {
    tmp <- which(unlist(lapply(lapply(x$bg, FUN = function(X) c(x$tree$edge[i, 1], x$tree$edge[i, 2]) %in% X), FUN = "all"), use.names = FALSE))
    if ((length(tmp) > 1)) {
      colo[i] <- tmp[which.min(unlist(lapply(x$bg, FUN = length)))]
    } else{
      colo[i] <- tmp
    }
  }
  ape::plot.phylo(x$tree, edge.color = colors[colo], ...)
}
#########Check bg list of nodes argument###########
checkbglist <- function(bg, tree){
  all_branches <- TRUE
  for (i in 1:length(tree$edge.length)) {
    if (!any(unlist(lapply(lapply(bg, FUN = function(j) tree$edge[i,] %in% j),FUN = "all")))) all_branches <- FALSE
  }
  return(all_branches)
}
patterns <- function(x) { 
  #x is the matrix of phyletic patterns...
  if (is.data.frame(x)) x <- as.matrix(x)
  wfn <- summary(as.factor(apply(x, 1, paste, collapse = " ")), maxsum = nrow(x))
  b <- attr(wfn, "names")
  ww <- unname(cbind(b, wfn))
  if (is.character(x[1, 1])) {
    databp_red1 <- na.omit(unlist(strsplit(ww[,1],split = " ")))
  } else {
    databp_red1 <- na.omit(as.numeric(unlist(strsplit(ww[,1],split = " "))))
  }
  databpr <- matrix(databp_red1, ncol = ncol(x), byrow = T)
  colnames(databpr) <- colnames(x)
  return(list(w = wfn, databp_red = databpr, names = colnames(x)))
}
estimaterates <- function(usertree = NULL, userphyl = NULL, matchtipstodata = FALSE, unobserved = NULL, alphabet = NULL, modelmat = NULL, bgtype = "listofnodes", bg = NULL, partition = NULL, ratevar = FALSE, nocat = 4, reversible = FALSE, numhessian = TRUE,  rootprob = NULL, rpvec = NULL, lowli = 0.001, upli = 100, ...) {
  ptm <- proc.time()
  if (ratevar != FALSE) ratevar <- match.arg(ratevar,c("discgamma","partitionspecificgamma"))
  bgtype <- match.arg(bgtype, c("listofnodes", "ancestornodes"))
  if (!is.null(rootprob)) rootprob <- match.arg(rootprob, c("equal", "stationary", "maxlik", "user"))
  if (is.character(modelmat)) modelmat <- match.arg(modelmat, c("ER","SYM","ARD","GTR","BD","BDER","BDSYM","BDARD"))
  if (is.null(usertree)) 
    stop("usertree option is required.")
  if (is.null(userphyl) | !(is.matrix(userphyl) | is.data.frame(userphyl))) 
    stop("userphyl option is required.")
  if (!is.null(unobserved) & !is.matrix(userphyl)) 
    stop("unobserved should be a matrix.")
  if (bgtype == "listofnodes" & !is.null(bg) & !is.list(bg))
    stop("bg should be a vector with listofnodes argument.")
  if (bgtype == "listofnodes" & !is.null(bg) & is.list(bg)) {
    if (!checkbglist(bg, usertree)) stop("Some branch is missing in bg.")
  }
  if (bgtype == "ancestornodes" & !is.null(bg) & !is.vector(bg))
    stop("bg should be a list with ancestornodes argument.")
  if (is.null(alphabet)) 
    stop("alphabet option is required.")
  if (is.null(modelmat)) 
    stop("modelmat option is required.")
  if (is.null(rootprob)) 
    stop("rootprob option is required.")
  if (!is.null(partition) & !is.list(partition)) 
    stop("partition option must be a list.")
  if (ratevar == "partitionspecificgamma" & is.null(partition)) 
    stop("partition option must be specified and must be a list.")
  if (!is.character(modelmat)) {
    if (nrow(modelmat) != ncol(modelmat))
      stop("A square matrix should be provided for the modelmat option.")
    if (nrow(modelmat) != length(alphabet))
      stop("Length of alphabet supplied does not equal the number of states in the supplied rate matrix.")
  }
  if (length(table(as.matrix(userphyl))) != length(alphabet))
    stop("Number of discrete states in the data provided does not equal the state space specified in the alphabet option.")
  if (nrow(userphyl) < ncol(userphyl))
    warning("More taxa than patterns detected. Note that rows should represent the discrete character patterns and columns the different taxa. If this format was used, you can safely ignore this message.")
  #############Libraries###################
  libload <- function(funcval) {
    if (requireNamespace(paste(funcval), quietly = TRUE) == FALSE) 
      stop("Please install package dependencies before continuing.")
    suppressMessages(loadNamespace(funcval))
  }
  libload("ape")
  libload("numDeriv")
  libload("Rcpp")
  #########Build substitution rate index matrix#####################
  buildindexmat <- function(type = modelmat, states = length(alphabet)) {
    if (type == "ER") {
      #equal rates matrix
      modelm <- matrix(1,states,states)
      diag(modelm) <- NA
    }
    if (type == "ARD") {
      #all rates different matrix
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      modelm[!is.na(modelm)] <- 1:(states*(states - 1))
    }
    if (type == "SYM") {
      #symmetric matrix
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      modelm[lower.tri(modelm)] <- 1:(states*(states - 1)/2)
      modelm <- modelm + t(modelm)
    }
    if (type == "GTR") {
      #symmetric matrix
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      modelm[lower.tri(modelm)] <- 1:(states*(states - 1)/2)
      modelm <- modelm + t(modelm)
    }
    if (type == "BD") {
      #birth death matrix
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      if (states > 2) {
        diag(modelm[-1, ]) <- 1
        diag(modelm[, -1]) <- 2
      } 
      if (states == 2) {
        modelm[2,1] <- 1
        modelm[1,2] <- 2
      }
    }
    if (type == "BDER") {
      #birth death equal rates matrix
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      if (states > 2) {
        diag(modelm[-1, ]) <- 1
        diag(modelm[, -1]) <- 1
      }
      if (states == 2) {
        modelm[2,1] <- modelm[2,1] <- 1
      }
    }
    if (type == "BDARD") {
      #birth death all rates different
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      if (states > 2) {
        diag(modelm[-1, ]) <- 1:(states - 1)
        diag(modelm[, -1]) <- (states):((states - 1)*2)
      }
      if (states == 2) {
        modelm[2, 1] <- 1
        modelm[1, 2] <- 2
      }
    }
    if (type == "BDSYM") {
      #birth death symmetric
      modelm <- matrix(0,states,states)
      diag(modelm) <- NA
      if (states > 2) {
        diag(modelm[-1, ]) <- 1:(states - 1)
        modelm <- modelm + t(modelm)
      } 
      if (states == 2) {
        modelm[2,1] <- modelm[2,1] <- 1
      }
    }
    return(modelm)
  }
  if (is.character(modelmat)) {
    if (modelmat == "GTR") {
      reversible <- TRUE
      rootprob <- "maxlik"
    }
  }
  if (is.character(modelmat)) modelmat <- buildindexmat()
  #########Transition rate matrix and substitution rate matrix#######
  TPM_taxa <- function(rates, ad, ti, rpin) {
    j <- unlist(lapply(lapply(bg, FUN = function(X) c(ad[1], ad[2]) %in% X), FUN = "all"), use.names = FALSE)
    if (sum(j) > 1) j[j][-which.min(unlist(lapply(bg, FUN = length)))] <- FALSE
    for (gf in 1:length(rates[, j])) {
      Q[modelmat == gf] <- rates[gf, j]
    }
    if (reversible & rootprob == "maxlik") Q <- sweep(Q, MARGIN = 2, FUN = "*", rpin, check.margin = FALSE)
    diag(Q) <- -.rowSums(Q, m = al, n = al, na.rm = TRUE)
    return(expm(Q * ti))
  }
  ###########Data Import#########
  if (!ape::is.binary.tree(usertree) | !ape::is.rooted(usertree)) {
    usertree <- ape::multi2di(usertree)
    cat("Tree either not binary or not rooted.", "\n", "Transformed using multi2di from packege ape.", "\n")
    cat("Labels of nodes of interest might change.\n")
  }
  tree1 <- usertree
  tree1 <- ape::reorder.phylo(tree1, order = "postorder")
  
  if (is.data.frame(userphyl)) userphyl <- as.matrix(userphyl)
  diag(modelmat) <- NA
  if (reversible & !isSymmetric(modelmat)) stop("modelmat should be symmetric to use with the reversible option.")
  if (reversible & rootprob != "maxlik") stop("Root probabilities must be estimated using likelihood to use with the reversible option. Use maxlik option for rootprob.")
  
  if (matchtipstodata) {
    fin <- userphyl[,pmatch(usertree$tip.label, colnames(userphyl))] 
    datab <- fin
  } else {
    datab <- userphyl
  }
  phyl <- nrow(datab)
  nooftaxa <- length(tree1$tip.label)
  al <- length(alphabet)
  Q <- modelmat
  if (bgtype == "listofnodes" & is.null(bg)) {
    bg <- list(c(1:(nooftaxa * 2 - 1)))
  } else if (bgtype == "listofnodes" & !is.null(bg)) {
    bg <- bg
  } else if (bgtype == "ancestornodes") {
    uncl <- sort(bg, decreasing = TRUE) #unique clades
    for (i in 1:length(uncl)) {
      if (!is.na(uncl[i])) {
        try <- phangorn::Ancestors(tree1, uncl[i], type = "all")
        res <- match(try, uncl)
        uncl[na.omit(res)] <- NA
      }
    }
    bg <- uncl <- na.omit(uncl)
    ########
    if (is.null(bg)) 
      stop("bg not specified.")
    bgo <- bg
    bg <- list()
    for (i in 1:length(bgo)) bg[[i]] <- c(bgo[i], phangorn::Descendants(tree1, bgo[i], type = "all"))
    bg[[length(bg) + 1]] <- c(setdiff(1:(nooftaxa * 2 - 1), c(unlist(bg, use.names = FALSE))), bgo)
  }
  ##############Parameters and Variables################
  nointnodes <- nooftaxa - 1
  tips <- 1:nooftaxa
  len_tips <- length(tips)
  databp <- datab
  csp <- length(bg) #how many different clade-specific parameters
  if (!is.null(partition)) {
    psp <- length(partition)
  } else {
    partition <- list(c(1:nrow(databp)))
    psp <- 1
  }
  if (ratevar == FALSE) 
    nocat <- 1
  w <- list()
  databp_red <- list()
  for (i in 1:psp) {
    temp <- patterns(x = databp[partition[[i]], ])
    w[[i]] <- temp$w
    databp_red[[i]] <- temp$databp_red
  }
  logll <- NULL
  results <- list()
  nodelist <- c(setdiff(tree1$edge[, 2], tips), len_tips + 1)
  F <- cbind(tree1$edge, tree1$edge.length)
  npar_mmat <- sort(unique(as.vector(modelmat[modelmat > 0])), na.last = NA)
  le_npar_mmat <- length(npar_mmat)
  rates <- array(data = NA, dim = c(le_npar_mmat, csp, psp))
  #partition specific parameters are in the third dimension with a matrix of 
  #clade-specific parameters in each matrix in the third dimension.
  ###############Log-likelihood function###############
  temp <- F[, 2]
  pweights <- function(comp = NULL, repar) {
    #repar are the new reparameterized variables
    #comp are the number of components for which the probability 
    #weights are being determined
    #The maximum likelihood estimates are found using the reparametrization
    #however, the standard errors are found numerically using the original
    #parameters.
    denom <- 1 + sum(unlist(lapply(X = repar, FUN = exp), use.names = FALSE))
    pwei <- c(unlist(lapply(X = repar, FUN = exp), use.names = FALSE), 1)/denom
    return(pwei)
  }
  quan <- seq(1/(2 * nocat), (2 * nocat - 1)/(2 * nocat), by = 2/(2 * nocat))
  le_csp <- le_npar_mmat * csp * psp
  if (ratevar == "partitionspecificgamma") {
    alpharates <- matrix(0, nrow = psp, ncol = nocat)
    rvp <- psp #rate var parameters
  } else {
    alpharates <- matrix(0, nrow = 1, ncol = nocat)
    rvp <- 1
  }
  bamsp <- function(previtval, model) { #branch and model specific assignment
    if (rootprob == "maxlik" & ratevar == FALSE) {
      rates[] <- previtval[1:(le_csp)]
      alpharates[1] <- 1
      return(list(rates = rates, alpharates = alpharates, iroot = previtval[(le_csp + 1):le_prev]))
    } else if (rootprob != "maxlik" & ratevar == "discgamma") {
      rates[] <- previtval[1:(le_csp)]
      alpha <- previtval[(le_csp + 1):le_prev]
      for (i in 1:rvp) {
        alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
          nocat
      }
      return(list(rates = rates, alpha = alpha, alpharates = alpharates))
    } else if (rootprob == "maxlik" & ratevar == "discgamma") {
      rates[] <- previtval[1:(le_csp)]
      alpha <- previtval[(le_csp + 1)]
      for (i in 1:rvp) {
        alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
          nocat
      }
      iroot <- previtval[(le_csp + 2):le_prev]
      return(list(rates = rates, alpha = alpha, iroot = iroot, alpharates = alpharates))
    } else if (rootprob != "maxlik" & ratevar == "partitionspecificgamma") {
      rates[] <- previtval[1:(le_csp)]
      alpha <- previtval[(le_csp + 1):(le_csp + psp)]
      for (i in 1:rvp) {
        alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
          nocat
      }
      return(list(rates = rates, alpha = alpha, alpharates = alpharates))
    } else if (rootprob == "maxlik" & ratevar == "partitionspecificgamma") {
      rates[] <- previtval[1:(le_csp)]
      alpha <- previtval[(le_csp + 1):(le_csp + psp)]
      for (i in 1:rvp) {
        alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
          nocat
      }
      iroot <- previtval[(le_csp + psp + 1):le_prev]
      return(list(rates = rates, alpha = alpha, iroot = iroot, alpharates = alpharates))
    } else {
      rates[] <- previtval
      alpharates[1] <- 1
      return(list(rates = rates, alpharates = alpharates))
    }
  }
  bamspforresult <- function(previtval, model, partype = NULL) { 
    #branch and model specific assignment
    if (rootprob == "maxlik" & ratevar == FALSE) {
      rates[] <- previtval[1:(le_csp)]
      if (partype == "par") {
        iroot <- pweights(comp = al, repar = (previtval[(le_csp + 1):le_prev]))
      } else {
        iroot <- previtval[(le_csp + 1):le_prev]
      }
      return(list(rates = rates, iroot = iroot))
    } else if (rootprob != "maxlik" & ratevar == "discgamma") {
      rates[] <- previtval[1:(le_csp)]
      if (partype == "se") {
        alpha <- previtval[(le_csp + 1):(le_csp + 1)]
        return(list(rates = rates, alpha = alpha))
      }
      if (partype == "par") {
        alpha <- previtval[(le_csp + 1):(le_csp + 1)]
        for (i in 1:rvp) {
          alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
            nocat
        }
        return(list(rates = rates, alpha = alpha, alpharates = alpharates))
      }
    } else if (rootprob == "maxlik" & ratevar == "discgamma") {
      rates[] <- previtval[1:(le_csp)]
      if (partype == "se") {
        alpha <- previtval[(le_csp + 1):(le_csp + 1)]
        iroot <- previtval[(le_csp + 1 + 1):le_prev]
        return(list(rates = rates, alpha = alpha, iroot = iroot))
      }
      if (partype == "par") {
        alpha <- previtval[(le_csp + 1):(le_csp + 1)]
        for (i in 1:rvp) {
          alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
            nocat
        }
        iroot <- pweights(comp = al, repar = (previtval[(le_csp + rvp + 1):le_prev]))
        return(list(rates = rates, alpha = alpha, alpharates = alpharates, iroot = iroot))
      }
      return(list(rates = rates, alpha = alpha, iroot = iroot))
    } else if (rootprob != "maxlik" & ratevar == "partitionspecificgamma") {
      rates[] <- previtval[1:(le_csp)]
      if (partype == "se") {
        alpha <- previtval[(le_csp + 1):(le_csp + psp)]
        return(list(rates = rates, alpha = alpha))
      }
      if (partype == "par") {
        alpha <- previtval[(le_csp + 1):(le_csp + psp)]
        for (i in 1:rvp) {
          alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
            nocat
        }
        return(list(rates = rates, alpha = alpha, alpharates = alpharates))
      }
    } else if (rootprob == "maxlik" & ratevar == "partitionspecificgamma") {
      rates[] <- previtval[1:(le_csp)]
      if (partype == "se") {
        alpha <- previtval[(le_csp + 1):(le_csp + psp)]
        iroot <- previtval[(le_csp + psp + 1):le_prev]
        return(list(rates = rates, alpha = alpha, iroot = iroot))
      }
      if (partype == "par") {
        alpha <- previtval[(le_csp + 1):(le_csp + psp)]
        for (i in 1:rvp) {
          alpharates[i, ] <- qgamma(quan, shape = alpha[i], scale = alpha[i])/sum(qgamma(quan, shape = alpha[i], scale = alpha[i])) * 
            nocat
        }
        iroot <- pweights(comp = al, repar = (previtval[(le_csp + psp + 1):le_prev]))
        return(list(rates = rates, alpha = alpha, alpharates = alpharates, iroot = iroot))
      }
      return(list(rates = rates, alpha = alpha, iroot = iroot))
    } else {
      rates[] <- previtval
      return(list(rates = rates))
    }
  }
  alseq <- 1:length(alphabet)
  Lixi <- matrix(0, nooftaxa + nointnodes, al)
  Lixi_init <- list()
  patlen <- unlist(lapply(X = w, FUN = length))
  for (i in 1:length(w)) {
    Lixi_init[[i]] <- rep(list(matrix(data = 0, nrow = nooftaxa + nointnodes, ncol = al)), length = patlen[i])
  } 
  for (i in 1:length(w)) {
    for (j in 1:patlen[i]) {
      for (u in alseq) {
        Lixi_init[[i]][[j]][which(databp_red[[i]][j, ] == alphabet[u]), u] <- 1
      }
    }
  }
  if (!is.null(unobserved)) {
    Lixi_corr_init <- list()
    patlen_corr <- nrow(unobserved)
    Lixi_corr_init <- rep(list(matrix(data = 0, nrow = nooftaxa + nointnodes, ncol = al)), length = patlen_corr)
      for (j in 1:patlen_corr) {
        for (u in alseq) {
          Lixi_corr_init[[j]][which(unobserved[j, ] == alphabet[u]), u] <- 1
        }
      }
  }
  tmp <- numeric(al)
  forav <- rep(list(modelmat), psp)
 
  ll <- function(model, pm_ll, rootp_ll, Lixi_in_ll) {
    tmp <- part_loopC(nocat, nodelist, al, tree1$edge[, 1], tree1$edge[, 2], pm_ll, Lixi_in_ll, len_tips)
    logll <- log(rowSums(sweep(tmp, MARGIN = 2, rootp_ll, `*`)))
    return(logll)
  }
  pm <- list()
  unlist_w = unlist(w)
  #################Objective function################
  totalll <- function(previtval, model, rtype = "estimates", ...) {
    ptbrun <- bamsp(previtval, model)
    ptb_loc <- ptbrun$rates
    rtype <- rtype
    rootp <- rootpfun(rootprob, rtype, currpar = ptbrun)
    if (ratevar == "partitionspecificgamma") {
      for (k in 1:psp) {
        if (le_npar_mmat == 1) {
          pm[[k]] <- lapply(X = 1:nocat, FUN = function(j) lapply(1:nrow(F), function(i) TPM_taxa(t(as.matrix(ptb_loc[, , k])), ad = F[i, 1:2], ti = F[i, 3] * ptbrun$alpharates[k, j], rootp[[k]])))
        } else {
          pm[[k]] <- lapply(X = 1:nocat, FUN = function(j) lapply(1:nrow(F), function(i) TPM_taxa(as.matrix(ptb_loc[, , k]), ad = F[i, 1:2], ti = F[i, 3] * ptbrun$alpharates[k, j], rootp[[k]])))
        }
      }
    } else {
      for (k in 1:psp) {
        if (le_npar_mmat == 1) {
          pm[[k]] <- lapply(X = 1:nocat, FUN = function(j) lapply(1:nrow(F), function(i) TPM_taxa(t(as.matrix(ptb_loc[, , k])), ad = F[i, 1:2], ti = F[i, 3] * ptbrun$alpharates[1, j], rootp[[k]])))
        } else {
          pm[[k]] <- lapply(X = 1:nocat, FUN = function(j) lapply(1:nrow(F), function(i) TPM_taxa(as.matrix(ptb_loc[, , k]), ad = F[i, 1:2], ti = F[i, 3] * ptbrun$alpharates[1, j], rootp[[k]])))
        }
      }
    }
    # pm is a list with the first component referring to the partitions, the second to the discrete gamma categories, and the third to the edge lengths.

    phy <- unlist_w * unlist(lapply(X = 1:psp, FUN = function(i) ll(model, pm[[i]], rootp = rootp[[i]], Lixi_in_ll = Lixi_init[[i]])), use.names = FALSE)
    if (!is.null(unobserved) & is.null(partition)) {
      corr <- unlist(lapply(X = 1:psp, FUN = function(i) ll(model, pm[[i]], rootp = rootp[[i]], Lixi_in_ll = Lixi_corr_init)), use.names = FALSE)
      return(-(sum(phy) - nrow(databp) * log(1 - sum(exp(corr))) ) )
    } else if (!is.null(unobserved) & !is.null(partition)) {
      corr <- lapply(X = 1:psp, FUN = function(i) ll(model, pm[[i]], rootp = rootp[[i]], Lixi_in_ll = Lixi_corr_init))
      return( -(sum(phy) - 
                 sum(unlist(lapply(X = partition, FUN = length), 
                                  use.names = FALSE) * 
                 unlist(lapply(X = 
                          lapply(X = 
                                   lapply(X = corr, FUN = exp), 
                                 FUN = sum), 
                        FUN = function(x) log(1 - x)) ) ) ) )
    } else{
    return(-(sum(phy)))
    }
    # When the likelihood is being conditioned on observable patterns only, and a partitioned analysis is being done, the multitude of lapplys are basically computing on lists the equivalent of log(1 - sum(exp(corr))) and then these are multiplied by a vector of number of observations in each partition.
    #-ve of the value because being minimized by default.
  }
  rootpfun <- function(rootprob, rtype, currpar) {
    if (rootprob == "maxlik") 
      irootprob <- currpar$iroot
    if (rootprob == "equal") {
      return(rep(list(rep(1/al, al)), psp))
    } else if (rootprob == "maxlik") {
      if (rtype == "estimates") {
        return(rep(list(pweights(comp = al, repar = irootprob)), psp))
      } else {
        # i.e. rtype == "se"
        return(rep(list(c(irootprob, (1 - sum(irootprob)))), psp))
      }
    } else if (rootprob == "stationary") {
      rates <- currpar$rates
      dimr <- dim(rates)
      ratesforav <- lapply(X = 1:dimr[3], FUN = function(x) .rowMeans(rates[, , x], m = dimr[1], n = dimr[2]))
      for (p in 1:psp) {
        for (gf in 1:length(ratesforav[[p]])) {
          forav[[p]][modelmat == gf] <- ratesforav[[p]][gf]
        }
      }
      for (p in 1:psp) diag(forav[[p]]) <- -.rowSums(forav[[p]], m = al, n = al, na.rm = TRUE)
      return(lapply(X = 1:psp, FUN = function(x) expm(forav[[x]] * 100)[1, ]))
      #            #Check whether all rows are equal to get proper equi freq
    } else if (rootprob == "user") {
      return(rep(list(rpvec), psp))
    }
  }
  indelinit <- 0.9
  ################Estimation################
  alphastart = 0.5
  options(digits = 7)
  if (rootprob == "maxlik" & ratevar == FALSE) {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat)), rep(0, al - 1)), df = 1 * csp * psp * le_npar_mmat + 
                                 al - 1, pb = 1, lower = c(rep(lowli, (csp * psp * le_npar_mmat)), rep(-100, al - 1)), upper = c(rep(upli, (csp * psp * le_npar_mmat)), 
                                                                                                                                 rep(100, al - 1)), model = "wop"))
  } else if (rootprob == "maxlik" & ratevar == "discgamma") {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat)), alphastart, rep(0, al - 1)), df = 1 * csp * psp * 
                                 le_npar_mmat + 1 + al - 1, pb = 1, lower = c(rep(lowli, (csp * psp * le_npar_mmat)), 0.01, rep(-100, al - 1)), upper = c(rep(upli, 
                                                                                                                                                              (csp * psp * le_npar_mmat)), 100, rep(100, al - 1)), model = "wop"))
  } else if (rootprob != "maxlik" & ratevar == "discgamma") {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat)), alphastart), df = 1 * csp * psp * le_npar_mmat + 
                                 1, pb = 1, lower = c(rep(lowli, (csp * psp * le_npar_mmat)), 0.01), upper = c(rep(upli, (csp * psp * le_npar_mmat)), 100), model = "wop"))
  } else if (rootprob == "maxlik" & ratevar == "partitionspecificgamma") {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat)), rep(alphastart, psp), rep(0, al - 1)), df = 1 * csp * 
                                 psp * le_npar_mmat + psp + al - 1, pb = 1, lower = c(rep(lowli, (csp * psp * le_npar_mmat)), rep(0.01, psp), rep(-100, al - 1)), 
                               upper = c(rep(upli, (csp * psp * le_npar_mmat)), rep(100, psp), rep(100, al - 1)), model = "wop"))
  } else if (rootprob != "maxlik" & ratevar == "partitionspecificgamma") {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat)), rep(alphastart, psp)), df = 1 * csp * psp * le_npar_mmat + 
                                 psp, pb = 1, lower = c(rep(lowli, (csp * psp * le_npar_mmat)), rep(0.01, psp)), upper = c(rep(upli, (csp * psp * le_npar_mmat)), 
                                                                                                                           rep(100, psp)), model = "wop"))
  } else {
    modelop <- list(wop = list(start = c(rep(indelinit, (csp * psp * le_npar_mmat))), df = 1 * csp * psp * le_npar_mmat, pb = 1, lower = rep(lowli, 
                                                                                                                                             (csp * psp * le_npar_mmat)), upper = rep(upli, (csp * psp * le_npar_mmat)), model = "wop"))
  }
  
  parvec <- NULL
  sevec <- NULL
  convcheck <- NULL
  modelnames <- "wop"
  modelmat2 <- NA
  for (i in modelnames) {
    le_prev <- length(modelop[[i]]$start)
    res <- nlminb(start = modelop[[i]]$start, objective = totalll, model = i, lower = modelop[[i]]$lower, upper = modelop[[i]]$upper, 
                  rtype = "estimates", ...)
    if (rootprob == "maxlik") {
      tempres <- res$par
      tempres[(length(tempres) - al + 2):length(tempres)] <- pweights(comp = al, repar = res$par[(length(res$par) - al + 2):length(res$par)])[1:(al - 1)]
      if (numhessian) 
        suppressWarnings(res$hessian <- numDeriv::hessian(func = totalll, model = i, x = tempres, rtype = "se"))
    } else {
      if (numhessian) 
        res$hessian <- numDeriv::hessian(func = totalll, model = i, x = res$par)
    }
    res$parsep <- bamspforresult(res$par, i, partype = "par")
    res$revcombined$rates <- array(data = NA, dim = c(sum(!is.na(modelmat) & modelmat != 0), csp, psp))
    if (reversible) {
      Q_out <- modelmat
      modelmat2 <- modelmat
      modelmat2[modelmat2 != 0 & !is.na(modelmat2)] <- 1:(sum(!is.na(modelmat) & modelmat != 0))
      for (pa in 1:psp) {
        for (ca in 1:csp) {
          for (gf in 1:length(res$parsep$rates[,ca,pa])) {
            Q_out[modelmat == gf] <- res$parsep$rates[,ca,pa][gf]
          }
              Q_out <- sweep(Q_out, MARGIN = 2, FUN = "*", res$parsep$iroot, check.margin = FALSE) 
              #iroot is never a list because only ever used with "maxlik", never "stationary". 
              #the sweep method is ideal. However, because we are only indexing !is.na items, the multiplication by the diagonal matrix comprising the root probabilities method might also work...
          res$revcombined$rates[,ca,pa] <- Q_out[which(!is.na(modelmat2) & modelmat2 != 0)]
        }
      }
      res$modelmat2 <- modelmat2
    }
    res$df <- modelop[[i]]$df
    convcheck <- c(convcheck, res$convergence)
    if (res$convergence == 0) {
      parvec <- c(parvec, res$par[1:le_csp]) #previtval[1:(le_csp)]
      if (numhessian) {
        res$se <- try(sqrt(diag(solve(res$hessian))))
        res$parsep$se <- try(bamspforresult(res$se, i, partype = "se"))
        sevec <- try(c(sevec, res$se))
      }
      llhat <- -round(res$objective, 3)
      res$AIC <- 2 * llhat - 2 * res$df
      res$BIC <- 2 * llhat - log(phyl) * res$df
    }
    results[[i]] <- res
  }
  if (lowli %in% parvec || upli %in% parvec) {
    cat("Estimated parameters on interval bounds.")
  }
  if (!all(is.finite(sevec)) & numhessian) {
    cat("Something is not right with the standard errors.")
    cat("Check Hessian matrix estimate.\n")
    cat("Consider calculating bootstrap errors (make sure to use numhessian=FALSE).\n")
  }
  
  timetaken <- proc.time() - ptm
  val <- list(call = match.call(), conv = convcheck, time = timetaken, bgtype = bgtype, bg = bg, results = results, tree = tree1, data_red = databp_red, alphabet = alphabet, reversible = reversible,  w = w, taxa = nooftaxa, phyl = phyl, rootprob = rootprob, modelmat = modelmat, modelnames = modelnames)
  class(val) <- "markophylo"
  return(invisible(val))
}