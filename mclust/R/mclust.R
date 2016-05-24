Mclust <- function(data, G = NULL, modelNames = NULL, prior = NULL, 
                   control = emControl(), initialization = NULL, 
                   warn = mclust.options("warn"), x = NULL, ...) 
{
  call <- match.call()
  data <- data.matrix(data)
  if(!is.null(x))
    if(!inherits(x, "mclustBIC"))
      stop("If provided, argument x must be an object of class 'mclustBIC'.")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name("mclustBIC")
  mc[[2]] <- data
  Bic <- eval(mc, parent.frame())
  G <- attr(Bic, "G")
  modelNames <- attr(Bic, "modelNames")
  Sumry <- summary(Bic, data, G = G, modelNames = modelNames)
  if(length(Sumry)==0) return()
  if(!(length(G) == 1)) 
    { bestG <- length(tabulate(Sumry$cl))
      if(warn) 
        { if(bestG == max(G) & warn) 
             warning("optimal number of clusters occurs at max choice")
          else if(bestG == min(G) & warn) 
               warning("optimal number of clusters occurs at min choice")
        }
  }
  oldClass(Sumry) <- NULL
  
  Sumry$bic <- Sumry$bic[1]
  Sumry$hypvol <- if(is.null(attr(Bic, "Vinv"))) 
                    as.double(NA) else 1/attr(Bic, "Vinv")
  # df <- (2*Sumry$loglik - Sumry$bic)/log(Sumry$n)
  df <- if(is.null(Sumry$modelName)) NULL 
        else with(Sumry, nMclustParams(modelName, d, G, 
                                       noise = (!is.na(hypvol)),
                                       equalPro = attr(Sumry, "control")$equalPro))

  ans <- c(list(call = call, data = data, BIC = Bic, df = df), Sumry)
  orderedNames <- c("call", "data", "modelName", 
                    "n", "d", "G", 
                    "BIC", "bic", "loglik", "df", 
                    "hypvol", "parameters", "z", 
                    "classification", "uncertainty")
  structure(ans[orderedNames], class = "Mclust")  
}

print.Mclust <- function(x, digits = getOption("digits"), ...)
{
  cat("\'", class(x)[1], "\' model object:\n", sep = "")
  G <- x$G
  if(G == 0 & !is.null(attr(x$BIC, "Vinv")))
  { cat(" best model: single noise component\n")
    return(invisible()) }
  M <- mclustModelNames(x$modelName)$type
  cat(" best model: ", M, " (", x$model, ") with ", G, " components\n", sep = "")
  invisible()
}

summary.Mclust <- function(object, parameters = FALSE, classification = FALSE, ...)
{
  # collect info
  G  <- object$G
  noise <- if(is.na(object$hypvol)) FALSE else object$hypvol
  pro <- object$parameters$pro
  if(is.null(pro)) pro <- 1
  names(pro) <- if(noise) c(seq_len(G),0) else seq(G)
  mean <- object$parameters$mean
  if(object$d > 1)
    { sigma <- object$parameters$variance$sigma }
  else
    { sigma <- rep(object$parameters$variance$sigmasq, object$G)[1:object$G]
      names(sigma) <- names(mean) }
  if(is.null(object$density))
    title <- paste("Gaussian finite mixture model fitted by EM algorithm")
  else
    title <- paste("Density estimation via Gaussian finite mixture modeling")
  #
  obj <- list(title = title, n = object$n, d = object$d, 
              G = G, modelName = object$modelName, 
              loglik = object$loglik, df = object$df, 
              bic = object$bic, icl = icl(object),
              pro = pro, mean = mean, variance = sigma,
              noise = noise,
              prior = attr(object$BIC, "prior"), 
              classification = object$classification, 
              printParameters = parameters, 
              printClassification = classification)
  class(obj) <- "summary.Mclust"
  return(obj)
}

print.summary.Mclust <- function(x, digits = getOption("digits"), ...)
{
  cat(rep("-", nchar(x$title)),"\n",sep="")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)),"\n",sep="")
  #
  if(is.null(x$modelName))
    { cat("\nMclust model with only a noise component") }
  else
    { cat("\nMclust ", x$modelName, " (", 
          mclustModelNames(x$modelName)$type, ") model with ", 
          x$G, ifelse(x$G > 1, " components", " component"),
          if(x$noise) "\nand a noise term", ":\n\n",
          sep = "") 
    }
  #
  if(!is.null(x$prior))
    { cat("Prior: ")
      cat(x$prior$functionName, "(", 
          paste(names(x$prior[-1]), x$prior[-1], sep = " = ", 
                collapse = ", "), ")", sep = "")
      cat("\n\n")
  }
  #
  tab <- data.frame("log-likelihood" = x$loglik, "n" = x$n, 
                    "df" = x$df, "BIC" = x$bic, "ICL" = x$icl, 
                    row.names = "")
  print(tab, digits = digits)
  #
  cat("\nClustering table:")
  print(table(factor(x$classification, 
                     levels = { l <- seq(x$G)
                                if(is.numeric(x$noise)) l <- c(l,0) 
                                l })),
        digits = digits)
  #
  if(x$printParameters)
  { cat("\nMixing probabilities:\n")
    print(x$pro, digits = digits)
    cat("\nMeans:\n")
    print(x$mean, digits = digits)
    cat("\nVariances:\n")
    if(x$d > 1) 
      { for(g in 1:x$G)
           { cat("[,,", g, "]\n", sep = "")
             print(x$variance[,,g], digits = digits) }
      }
    else print(x$variance, digits = digits)
    if(x$noise)
      { cat("\nHypervolume of noise component:\n")
        cat(signif(x$noise, digits = digits), "\n") }
  }
  if(x$printClassification)
  { cat("\nClassification:\n")
    print(x$classification, digits = digits)
  }
  #
  invisible(x)
}

# old version, wrong df, data not needed, inefficient because compute dens on original scale
# logLik.Mclust <- function(object, data, ...)
# {
#   if(!missing(data))
#     object$data <- data.matrix(data)
#   par <- object$parameters
#   df <- with(object, (G-1) + G*d + nVarParams(modelName, d = d, G = G))
# 
#   #  l <- matrix(as.double(NA), object$n, object$G)
#   #  for(k in seq(object$G))
#   #     { l[,k] <- par$pro[k] * dmvnorm(data, par$mean[,k],
#   #                                           par$variance$sigma[,,k]) }
#   #  l <- sum(log(rowSums(l)))
#   
#   l <- sum(log(do.call("dens", object)))
#   attr(l, "nobs") <- object$n
#   attr(l, "df") <- df
#   class(l) <- "logLik"
#   return(l)
# }

logLik.Mclust <- function(object, ...)
{
  par <- object$parameters
  # df <- with(object, (G-1) + G*d + nVarParams(modelName, d = d, G = G))
  df <- with(object, nMclustParams(modelName, d, G, 
                                   noise = (!is.na(hypvol)), 
                                   equalPro = attr(BIC, "control")$equalPro))
  l <- sum(do.call("dens", c(object, logarithm = TRUE)))
  attr(l, "nobs") <- object$n
  attr(l, "df") <- df
  class(l) <- "logLik"
  return(l)
}

predict.Mclust <- function(object, newdata, ...)
{
  if(!inherits(object, "Mclust")) 
    stop("object not of class \"Mclust\"")
  if(missing(newdata))
    { newdata <- object$data }
  prior <- object$parameters$pro
  object$data <- newdata
  # old
  # z <- do.call("cdens", object)
  # z <- sweep(z, MARGIN = 1, FUN = "/", STATS = apply(z, 1, max))
  # z <- sweep(z, MARGIN = 2, FUN = "*", STATS = prior/sum(prior))
  # z <- sweep(z, MARGIN = 1, STATS = apply(z, 1, sum), FUN = "/")
  # new: more efficient and accurate
  z <- do.call("cdens", c(object, list(logarithm = TRUE)))
  z <- sweep(z, MARGIN = 2, FUN = "+", STATS = log(prior/sum(prior)))
  z <- sweep(z, MARGIN = 1, FUN = "-", STATS = apply(z, 1, logsumexp))
  z <- exp(z)
  #
  cl <- apply(z, 1, which.max)
  out <- list(classification = cl, z = z)
  return(out) 
}

#############################################################################

EMclust <- function(data, G = NULL, modelNames = NULL, prior = NULL, control = emControl(), initialization = list(hcPairs=NULL, subset=NULL, noise=NULL), Vinv = NULL, warn = FALSE, x = NULL, ...)
{
  if (!is.null(x)) {
    if (!missing(prior) || !missing(control) || 
          !missing(initialization) || !missing(Vinv)) 
      stop("only G and modelNames may be specified as arguments when x is supplied")
    prior <- attr(x,"prior") 
    control <- attr(x,"control")
    initialization <- attr(x,"initialization")
    Vinv <- attr(x,"Vinv")
    warn <- attr(x,"warn")
  }
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  if (is.null(x)) {
    if (is.null(modelNames)) {
      if (d == 1) {
        modelNames <- c("E", "V")
      }
      else {
        modelNames <- mclust.options("emModelNames")
        if (n <= d) {          
          # select only spherical and diagonal models
          m <- match(modelNames, c("EII", "VII", "EEI", "VEI", "EVI", "VVI"),
                     nomatch = 0)
          modelNames <- modelNames[m]
        }
      }
    }
    if (is.null(G)) {
      G <- if (is.null(initialization$noise)) 1:9 else 0:9
    }
    else {
      G <- sort(as.numeric(G))
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if (is.null(G)) G <- Glabels
    if (is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if (all(Gmatch) && all(Mmatch)) {
      attr( x, "G") <- as.numeric(G)
      attr( x, "modelNames") <- modelNames
      attr( x, "returnCodes") <- 
        attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
      return(x[as.character(G),modelNames,drop=FALSE])
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  if (any(as.logical(as.numeric(G))) < 0) {
    if (is.null(initialization$noise)) {
      stop("G must be positive")
    }
    else {
      stop("G must be nonnegative")
    }
  }
  if (d == 1 && any(nchar(modelNames) > 1)) {
    Emodel <- any(sapply(modelNames, function(x)
      charmatch("E", x, nomatch = 0)[1]) == 1)
    Vmodel <- any(sapply(modelNames, function(x)
      charmatch("V", x, nomatch = 0)[1]) == 1)
    modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if (!is.null(x)) {
    BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
    RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
    BIC <- BIC[as.character(G),modelNames,drop=FALSE]
    RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1",] == EMPTY]) {
        out <- mvn(modelName = mdl, data = data, prior = prior)
        BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                             n = n, d = d, G = 1, equalPro = FALSE)
        RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
      #######################################################
      # all data in initial hierarchical clustering phase
      #######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], data = data)
          }
          else {
            hcPairs <- hc(modelName = "EII", data = data)
          } 
        }
        else {
          hcPairs <- NULL 
          #   hcPairs <- hc(modelName = "E", data = data)
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          z <- unmap(clss[, g]) 
        }
        else {
          z <- unmap(qclass( data, as.numeric(g)))
        }
        for (modelName in modelNames[BIC[g,] == EMPTY]) {
          out <- me(modelName = modelName, data = data, z = z, 
                    prior = prior, control = control, warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik,
                                   n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
        }
      }
    }
    else {
      ######################################################
      # initial hierarchical clustering phase on a subset
      ######################################################
      if (is.logical(initialization$subset)) 
        initialization$subset <- (1:n)[initialization$subset]
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(data = data[initialization$subset,  ],
                          modelName = mclust.options("hcModelNames")[1])
          }
          else {
            hcPairs <- hc(data = data[initialization$subset,],
                          modelName = "EII")
          }
        }
        else {
          hcPairs <- NULL
          #    hcPairs <- hc(data = data[initialization$subset],
          #                  modelName = "E")
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          z <- unmap(clss[, g]) 
        }
        else {
          z <- unmap(qclass(data[initialization$subset], as.numeric(g)))
        }
        dimnames(z) <- list(as.character(initialization$subset), NULL)
        for (modelName in modelNames[!is.na(BIC[g,])]) {
          ms <- mstep(modelName = modelName, z = z, 
                      data = as.matrix(data)[initialization$subset,  ],
                      prior = prior, control = control, warn = warn)
          #
          #  ctrl <- control
          #  ctrl$itmax[1] <- 1
          #  ms <- me(modelName = modelName, data = as.matrix(data)[
          #           initialization$subset,  ], z = z, prior = prior, control = ctrl)
          #
          es <- do.call("estep", c(list(data = data, warn = warn), ms))
          out <- me(modelName = modelName, data = data, z = es$z, 
                    prior = prior, control = control, warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik,
                                   n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
        }
      }
    }
  }
  else {
    ######################################################
    # noise case
    ######################################################
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    noise <- initialization$noise
    if (!is.logical(noise))
      noise <- as.logical(match(1:n, noise, nomatch = 0))
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
                         initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(data = data[!noise,],
                        modelName = mclust.options("hcModelNames")[1])
        }
        else {
          hcPairs <- hc(data = data[!noise,], modelName = "EII")
        }
      }
      else {
        hcPairs <- NULL 
        #    hcPairs <- hc(data = data[!noise], modelName = "E")
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
      z[] <- 0
      k <- as.numeric(g)
      if (d > 1 || !is.null(hcPairs)) {
        z[!noise, 1:k] <- unmap(clss[, g])
      }
      else {
        z[!noise, 1:k] <- unmap(qclass(data[!noise]))
      }
      z[noise, k+1] <- 1
      K <- 1:(k+1) 
      for (modelName in modelNames[BIC[g,] == EMPTY]) {
        out <- me(modelName = modelName, data = data, z = z[, K], 
                  prior = prior, control = control, Vinv = Vinv, warn = warn)
        BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                 n = n, d = d, G = k, 
                                 noise = TRUE, 
                                 equalPro = control$equalPro)
        RET[g, modelName] <- attr(out, "returnCode")
      }
    }
  }
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            returnCodes = RET, class = "mclustBIC")
}

# EMclust <- function(...) .Defunct("mclustBIC", PACKAGE = "mclust")
#############################################################################

mclustBIC <- function(data, G = NULL, modelNames = NULL, 
                      prior = NULL, control = emControl(), 
                      initialization = list(hcPairs=NULL, subset=NULL, noise=NULL),  
                      Vinv = NULL, warn = mclust.options("warn"), 
                      x = NULL, ...)
{
  if(!is.null(x)) 
    { if(!missing(prior) || !missing(control) || 
         !missing(initialization) || !missing(Vinv))
         stop("only G and modelNames may be specified as arguments when x is supplied")
      prior <- attr(x,"prior") 
      control <- attr(x,"control")
      initialization <- attr(x,"initialization")
      Vinv <- attr(x,"Vinv")
      warn <- attr(x,"warn")
  }
  dimData <- dim(data)
  oneD <- (is.null(dimData) || length(dimData[dimData > 1]) == 1)
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) 
    { data <- drop(as.matrix(data))
      n <- length(data)
      d <- 1
  }
  else {
      data <- as.matrix(data)
      n <- nrow(data)
      d <- ncol(data)
  }
  if(is.null(x)) 
    { if(is.null(modelNames)) 
        { if(d == 1) 
            { modelNames <- c("E", "V") }
          else {
              modelNames <- mclust.options("emModelNames")
              if(n <= d) 
                { # select only spherical and diagonal models
                  m <- match(modelNames, c("EII", "VII", "EEI", 
                                           "VEI", "EVI", "VVI"),
                             nomatch = 0)
                  modelNames <- modelNames[m]
                }
          }
    }
    if(!is.null(prior))
      { # remove models not available with prior
        modelNames <- setdiff(modelNames, c("EVE","VEE","VVE","EVV"))
    }
    if(is.null(G)) 
      { G <- if (is.null(initialization$noise)) 1:9 else 0:9 }
    else {
        G <- sort(as.integer(unique(G)))
    }
    if(is.null(initialization$noise)) 
      { if (any(G > n)) G <- G[G <= n] }
    else {
        noise <- initialization$noise
        if(!is.logical(noise)) 
          { if(any(match(noise, 1:n, nomatch = 0) == 0))
               stop("numeric noise must correspond to row indexes of data")
            noise <- as.logical(match(1:n, noise, nomatch = 0))
        }
        initialization$noise <- noise
        nnoise <- sum(as.numeric(noise))
        if(any(G > (n-nnoise))) G <- G[G <= n-nnoise]
    }
    if(!is.null(initialization$subset)) 
      { subset <- initialization$subset
        if(is.logical(subset)) subset <- which(subset)
        n <- length(subset)
        if(any(G > n)) G <- G[G <= n]
    }
    Gall <- G
    Mall <- modelNames
  }
  else {
    Glabels <- dimnames(x)[[1]]
    Mlabels <- dimnames(x)[[2]]
    if(is.null(G)) G <- Glabels
    if(is.null(modelNames)) modelNames <- Mlabels
    Gmatch <- match(as.character(G), Glabels, nomatch = 0)
    Mmatch <- match(modelNames, Mlabels, nomatch = 0)
    if(all(Gmatch) && all(Mmatch)) 
      { out <- x[as.character(G),modelNames,drop=FALSE]
        mostattributes(out) <- attributes(x)
        attr(out, "dim") <- c(length(G), length(modelNames))
        attr(out, "dimnames") <- list(G, modelNames)
        attr(out, "G") <- as.numeric(G)
        attr(out, "modelNames") <- modelNames
        attr(out, "returnCodes") <- 
          attr(x, "returnCodes")[as.character(G),modelNames,drop=FALSE]
        return(out)
    }
    Gall <- sort(as.numeric(unique(c(as.character(G), Glabels))))
    Mall <- unique(c(modelNames, Mlabels))
  }
  
  if(any(as.logical(as.numeric(G))) < 0) 
    { if(is.null(initialization$noise)) 
        { stop("G must be positive") }
    else {
          stop("G must be nonnegative")
    }
  }
  
  if(d == 1 && any(nchar(modelNames) > 1)) 
    { Emodel <- any(sapply(modelNames, function(x)
                    charmatch("E", x, nomatch = 0)[1]) == 1)
      Vmodel <- any(sapply(modelNames, function(x)
                    charmatch("V", x, nomatch = 0)[1]) == 1)
      modelNames <- c("E", "V")[c(Emodel, Vmodel)]
  }
  l <- length(Gall)
  m <- length(Mall)
  EMPTY <- -.Machine$double.xmax
  BIC <- RET <- matrix(EMPTY, nrow = l, ncol = m, 
                       dimnames = list(as.character(Gall), as.character(Mall)))
  if(!is.null(x)) 
    { BIC[dimnames(x)[[1]],dimnames(x)[[2]]] <- x
      RET[dimnames(x)[[1]],dimnames(x)[[2]]] <- attr(x, "returnCodes")
      BIC <- BIC[as.character(G),modelNames,drop=FALSE]
      RET <- RET[as.character(G),modelNames,drop=FALSE]
  }
  G <- as.numeric(G)
  Glabels <- as.character(G)
  Gout <- G
  if (is.null(initialization$noise)) {
    if (G[1] == 1) {
      for (mdl in modelNames[BIC["1",] == EMPTY]) {
        out <- mvn(modelName = mdl, data = data, prior = prior)
        BIC["1", mdl] <- bic(modelName = mdl, loglik = out$loglik, 
                             n = n, d = d, G = 1, equalPro = FALSE)
        RET["1", mdl] <- attr(out, "returnCode")
      }
      if (l == 1) {
        BIC[BIC == EMPTY] <- NA
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, initialization = initialization, 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$subset)) {
      #######################################################
      # all data in initial hierarchical clustering phase
      #######################################################
      if (is.null(initialization$hcPairs)) { 
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(data = data, 
                          modelName = mclust.options("hcModelNames")[1])
          }
          else {
            hcPairs <- hc(data = data, modelName = "EII")
          } 
        }
        else {
          hcPairs <- NULL 
          #   hcPairs <- hc(data = data, modelName = "E")
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs))  clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[,g]
        }
        else {
          cl <- qclass( data, as.numeric(g))
        }
        z <- unmap(cl, groups = 1:max(cl))
        if(any(apply( z, 2, max) == 0) & warn) 
          { #  missing groups
            if(warn) warning("there are missing groups")
            small <- sqrt(.Machine$double.neg.eps)
            z[z < small] <- small
            z <-  t(apply( z, 1, function(x) x/sum(x)))
        }
#       for(modelName in modelNames[BIC[g,] == EMPTY]) {
#       LS:
        for(modelName in na.omit(modelNames[BIC[g,] == EMPTY])) {  
          out <- me(modelName = modelName, data = data, z = z, 
                    prior = prior, control = control, warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik,
                                   n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
        }
      }
    }
    else {
      subset <- initialization$subset
      if (is.logical(subset)) subset <- which(subset)
      ######################################################
      # initial hierarchical clustering phase on a subset
      ######################################################
      if (is.null(initialization$hcPairs)) {
        if (d != 1) {
          if (n > d) {
            hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                          data = data[subset,])
          }
          else {
            hcPairs <- hc(modelName = "EII", 
                          data = data[subset,])
          }
        }
        else {
          hcPairs <- NULL
          #    hcPairs <- hc(modelName = "E", 
          #                  data = data[subset])
        }
      }
      else hcPairs <- initialization$hcPairs
      if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
      for (g in Glabels) {
        if (d > 1 || !is.null(hcPairs)) {
          cl <- clss[, g]
        }
        else {
          cl <- qclass(data[subset], as.numeric(g))
        }
        z <- unmap(cl, groups = 1:max(cl))
        if(any(apply( z, 2, max) == 0) & warn) 
          { #  missing groups
            if(warn) warning("there are missing groups")
            small <- sqrt(.Machine$double.neg.eps)
            z[z < small] <- small
            z <-  t(apply( z, 1, function(x) x/sum(x)))
        }
        for (modelName in modelNames[!is.na(BIC[g,])]) {
          ms <- mstep(modelName = modelName, z = z, 
                      data = as.matrix(data)[initialization$subset,  ],
                      prior = prior, control = control, warn = warn)
          #
          #  ctrl <- control
          #  ctrl$itmax[1] <- 1
          #  ms <- me(modelName = modelName, data = as.matrix(data)[
          #           initialization$subset,  ], z = z, prior = prior, control = ctrl)
          #
          es <- do.call("estep", c(list(data = data, warn = warn), ms))
          out <- me(modelName = modelName, data = data, z = es$z, 
                    prior = prior, control = control, warn = warn)
          BIC[g, modelName] <- bic(modelName = modelName, 
                                   loglik = out$loglik,
                                   n = n, d = d, G = as.numeric(g), 
                                   equalPro = control$equalPro)
          RET[g, modelName] <- attr(out, "returnCode")
        }
      }
    }
  }
  else {
    ######################################################
    # noise case
    ######################################################
    noise <- initialization$noise
    if (!is.null(initialization$subset)) 
      stop("subset option not implemented with noise")
    if (is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
    if (!G[1]) {
      hood <- n * logb(Vinv)
      BIC["0",  ] <- 2 * hood - logb(n)
      if (l == 1) {
        return(structure(BIC, G = G, modelNames = modelNames, prior = prior, 
                         control = control, 
                         initialization = list(hcPairs = hcPairs, subset = initialization$subset), 
                         warn = warn, n = n, d = d, oneD = oneD,
                         returnCodes = RET, class =  "mclustBIC"))
      }
      G <- G[-1]
      Glabels <- Glabels[-1]
    }
    if (is.null(initialization$hcPairs)) {
      if (d != 1) {
        if (n > d) {
          hcPairs <- hc(modelName = mclust.options("hcModelNames")[1], 
                        data = data[!noise,  ])
        }
        else {
          hcPairs <- hc(modelName = "EII", data = data[!noise,  ])
        }
      }
      else {
        hcPairs <- NULL 
        #    hcPairs <- hc(modelName = "E", data = data[!noise])
      }
    }
    else hcPairs <- initialization$hcPairs
    if (d > 1 || !is.null(hcPairs)) clss <- hclass(hcPairs, G)
    z <- matrix(0, n, max(G) + 1)
    for (g in Glabels) {
      z[] <- 0
      k <- as.numeric(g)
      if (d > 1 || !is.null(hcPairs)) {
        cl <- clss[, g]
      }
      else {
        cl <- qclass(data[!noise], k = k)
      }
      z[!noise,1:k] <- unmap(cl, groups = 1:max(cl))
      if(any(apply( z[!noise,1:k,drop=FALSE], 2, max) == 0) & warn) 
        { #           missing groups
          if(warn) warning("there are missing groups")         
          z[!noise,1:k] <- max( z[!noise,1:k], sqrt(.Machine$double.neg.eps))
          z[!noise,1:k] <- apply( z[!noise,1:k,drop=FALSE], 1, function(z) z/sum(z))
      }
      z[noise, k+1] <- 1
      K <- 1:(k+1) 
#     for (modelName in modelNames[BIC[g,] == EMPTY]) {
#     LS:
      for (modelName in na.omit(modelNames[BIC[g,] == EMPTY])) {  
        out <- me(modelName = modelName, data = data, z = z[, K], 
                  prior = prior, control = control, Vinv = Vinv, warn = warn)
        BIC[g, modelName] <- bic(modelName = modelName, loglik = out$loglik, 
                                 n = n, d = d, G = k, 
                                 noise = TRUE, 
                                 equalPro = control$equalPro)
        RET[g, modelName] <- attr(out, "returnCode")
      }
    }
  }
  
  if(!is.null(prior) & any(is.na(BIC)))
    warning("The presence of BIC values equal to NA is likely due to one or more of the mixture proportions being estimated as zero, so that the model estimated reduces to one with a smaller number of components.")  
  
  structure(BIC, G = Gout, modelNames = modelNames, prior = prior, 
            control = control, 
            initialization = list(hcPairs = hcPairs, 
                                  subset = initialization$subset,
                                  noise = initialization$noise), 
            Vinv = Vinv, warn = warn, n = n, d = d, oneD = oneD,
            criterion = "BIC", returnCodes = RET, 
            class = "mclustBIC")
}

pickBIC <- function(x, k = 3)
{
  if(!is.matrix(x))
    { warning("sorry, the pickBIC function cannot be applied to the provided argument!")
    return() }
  Glabels <- dimnames(x)[[1]]
  modelNames <- dimnames(x)[[2]]
  mis <- is.na(x)
  if(all(mis)) 
    { warning("none of the selected models could be fitted")
      return(rep(NA,k))
  }
  x[mis] <-  - .Machine$double.xmax
  x <- data.frame(as.vector(x), Glabels[as.vector(row(x))], 
                  modelNames[as.vector(col(x))])
  # x <- x[rev(order(x[,1])),] 
  # order by including first simpler models if ties are present
  x <- x[order(-x[, 1], x[,2], x[,3]),]
  namesx <- apply(x[,-1,drop = FALSE], 1, function(z) 
    paste(as.character(z[2]), as.character(z[1]), sep = ","))
  k <- min(k, nrow(x))
  x <- x[1:k,1]
  x[x ==  - .Machine$double.xmax] <- NA
  namesx <- namesx[1:k]
  namesx[is.na(x)] <- " "
  names(x) <- namesx
  x
}

print.mclustBIC <- function(x, pick = 3, ...)
{
  subset <- !is.null(attr(x, "subset"))
  oldClass(x) <- attr(x, "args") <- NULL
  attr(x, "criterion") <- NULL
  attr(x, "control") <- attr(x, "initialization") <- NULL
  attr(x, "oneD") <- attr(x, "warn") <- attr(x, "Vinv") <- NULL
  attr(x, "prior") <- attr(x, "G") <- attr(x, "modelNames") <- NULL
  ret <- attr(x, "returnCodes") == -3
  n <- attr(x, "n")
  d <- attr(x, "d")
  attr(x, "returnCodes") <- attr(x, "n") <- attr(x, "d") <- NULL
  
  cat("Bayesian Information Criterion (BIC):\n")
  NextMethod("print")
  cat("\n")
  cat("Top", pick, "models based on the BIC criterion:\n")
  print(pickBIC(x, pick), ...)
  invisible()
}

mclustModel <- function(data, BICvalues, G=NULL, modelNames=NULL, ...)
{
  mc <- match.call(expand.dots = FALSE)
  if (is.null(attr(BICvalues,"initialization")$noise)) {
    mc[[1]] <- as.name("summaryMclustBIC")
  }
  else {
    mc[[1]] <- as.name("summaryMclustBICn")
  }
  nm <- names(mc)
  mc[1:3] <- mc[c(1,3,2)]
  nm[1:3] <- nm[c(1,3,2)]
  nm[nm == "BICvalues"] <- "object" 
  names(mc) <- nm
  ans <- eval(mc, parent.frame())
  ans$classification <- ans$uncertainty <- NULL
  attr( ans, "bestBICvalues") <- NULL
  attr( ans, "prior") <- NULL
  attr( ans, "control") <- NULL
  attr( ans, "initialization") <- NULL
  oldClass(ans) <- "mclustModel"
  ans
}

mclustModelNames <- function(model)
{
  type <- switch(EXPR = as.character(model),
                 "E" = "univariate, equal variance",
                 "V" = "univariate, unequal variance",
                 "EII" = "spherical, equal volume",
                 "VII" = "spherical, varying volume",
                 "EEI" = "diagonal, equal volume and shape",
                 "VEI" = "diagonal, equal shape",
                 "EVI" = "diagonal, equal volume, varying shape",
                 "VVI" = "diagonal, varying volume and shape",
                 "EEE" = "ellipsoidal, equal volume, shape and orientation",
                 "EVE" = "ellipsoidal, equal volume and orientation",
                 "VEE" = "ellipsoidal, equal shape and orientation",
                 "VVE" = "ellipsoidal, equal orientation",
                 "EEV" = "ellipsoidal, equal volume and shape",
                 "VEV" = "ellipsoidal, equal shape",
                 "EVV" = "ellipsoidal, equal volume",
                 "VVV" = "ellipsoidal, varying volume, shape, and orientation",
                 "X"   = "univariate normal",
                 "XII" = "spherical multivariate normal",
                 "XXI" = "diagonal multivariate normal",
                 "XXX" = "ellipsoidal multivariate normal",
                 warning("invalid model"))
  return(list(model = model, type = type))
}

summary.mclustBIC <- function(object, data, G, modelNames, ...)
{
  mc <- match.call(expand.dots = FALSE)
  if(missing(data)) 
    { if(!missing(G)) 
        object <- object[rownames(object) %in% G,,drop=FALSE]
      if(!missing(modelNames)) 
        object <- object[,colnames(object) %in% modelNames,drop=FALSE]
      ans <- pickBIC(object, ...)
      class(ans) <- "summary.mclustBIC" 
  } else
  { if(is.null(attr(object,"initialization")$noise)) 
      { mc[[1]] <- as.name("summaryMclustBIC") }
    else 
      { mc[[1]] <- as.name("summaryMclustBICn") }
    warn <- attr(object, "warn")
    ans <- eval(mc, parent.frame())
    Glabels <- dimnames(object)[[1]]
    if(length(Glabels) != 1 && (!missing(G) && length(G) > 1)) 
      { Grange <- range(as.numeric(Glabels))
        if(match(ans$G, Grange, nomatch = 0) & warn)
          warning("best model occurs at the min or max of number of components considered!")
    }
  }
  ans
}


summaryMclustBIC <- function (object, data, G = NULL, modelNames = NULL, ...) 
{
  dimData <- dim(data)
  oneD <- (is.null(dimData) || length(dimData[dimData > 1]) == 1)
  if (!oneD && length(dimData) != 2) 
    stop("data must be a vector or a matrix")
  if (oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  initialization <- attr(object, "initialization")
  hcPairs <- initialization$hcPairs
  subset <- initialization$subset
  prior <- attr(object, "prior")
  control <- attr(object, "control")
  warn <- attr(object, "warn")
  oldClass(object) <- NULL
  attr(object, "prior") <- attr(object, "warn") <- NULL
  attr(object, "modelNames") <- attr(object, "oneD") <- NULL
  attr(object, "initialization") <- attr(object, "control") <- NULL
  d <- if (is.null(dim(data))) 1 else ncol(data)
  if(is.null(G)) 
    G <- dimnames(object)[[1]]
  if(is.null(modelNames)) 
    modelNames <- dimnames(object)[[2]]
  bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
  if(all(is.na(bestBICs))) 
    { return(structure(NULL, bestBICvalues = bestBICs, prior = prior, 
                       control = control, initialization = initialization, 
                       class = "summary.mclustBIC")) 
  }
  temp <- unlist(strsplit(names(bestBICs)[1], ","))
  bestModel <- temp[1]
  G <- as.numeric(temp[2])
  if(G == 1) 
    { out <- mvn(modelName = bestModel, data = data, prior = prior)
      ans <- c(list(bic = bestBICs, 
                    z = unmap(rep(1,n)),
                    classification = rep(1, n), 
                    uncertainty = rep(0, n)), 
               out)
    orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", 
                      "parameters", "z", "classification", "uncertainty")
    return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                     prior = prior, control = control, 
                     initialization = initialization, 
                     class = "summary.mclustBIC"))
  }
  
  if(is.null(subset)) 
  {
    if(d > 1 || !is.null(hcPairs))
      { z <- unmap(hclass(hcPairs, G)) }
    else 
      { z <- unmap(qclass(data, G), groups = 1:G) }
    out <- me(modelName = bestModel, data = data, z = z, 
              prior = prior, control = control, warn = warn)
    if(sum((out$parameters$pro - colMeans(out$z))^2) > 
           sqrt(.Machine$double.eps))
      { # perform extra M-step and update parameters
        out$parameters <- mstep(modelName = bestModel, data = data, 
                                z = out$z, prior = prior, 
                                warn = warn)$parameters 
    }
  }
  else 
  {
    if(d > 1 || !is.null(hcPairs)) 
      { z <- unmap(hclass(hcPairs, G)) }
    else 
      { z <- unmap(qclass(data[subset], G)) }
    ms <- mstep(modelName = bestModel, prior = prior, z = z, 
                data = as.matrix(data)[subset,], control = control, 
                warn = warn)
    es <- do.call("estep", c(list(data = data), ms))
    out <- me(modelName = bestModel, data = data, z = es$z, 
              prior = prior, control = control, warn = warn)
    # perform extra M-step and update parameters
    out$parameters <- mstep(modelName = bestModel, data = data, 
                            z = out$z, prior = prior, 
                            warn = warn)$parameters
  }
  obsNames <- if (is.null(dim(data))) names(data) else dimnames(data)[[1]]
  classification <- map(out$z, warn = warn)
  uncertainty <- 1 - apply(out$z, 1, max)
  names(classification) <- names(uncertainty) <- obsNames
  
  ans <- c(list(bic = bic(bestModel, out$loglik, out$n, out$d, out$G, 
                          noise = FALSE, equalPro = control$equalPro),
                # bic = as.vector(bestBICs[1]), 
                classification = classification, 
                uncertainty = uncertainty), 
           out)
  orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", 
                    "parameters", "z", "classification", "uncertainty")
  structure(ans[orderedNames], bestBICvalues = bestBICs, prior = prior, 
            control = control, initialization = initialization, 
            class = "summary.mclustBIC")
}

summaryMclustBICn <- function(object, data, G = NULL, modelNames = NULL, ...)
{
  dimData <- dim(data)
  oneD <- is.null(dimData) || length(dimData[dimData > 1]) == 1
  if(!oneD && length(dimData) != 2)
    stop("data must be a vector or a matrix")
  if(oneD) {
    data <- drop(as.matrix(data))
    n <- length(data)
    d <- 1
  }
  else {
    data <- as.matrix(data)
    n <- nrow(data)
    d <- ncol(data)
  }
  initialization <- attr(object, "initialization")
  hcPairs <- initialization$hcPairs
  noise <- initialization$noise
  if(!is.logical(noise))
    noise <- as.logical(match(1:n, noise, nomatch = 0))
  prior <- attr(object, "prior")
  control <- attr(object, "control")
  warn <- attr(object, "warn")
  Vinv <- attr(object, "Vinv")
  oldClass(object) <- NULL
  attr(object, "control") <- attr(object, "initialization") <- NULL
  attr(object, "prior") <- attr(object, "Vinv") <- NULL
  attr(object, "warn") <- NULL
  ##
  if (is.null(G)) G <- dimnames(object)[[1]]
  if (is.null(modelNames)) modelNames <- dimnames(object)[[2]]
  bestBICs <- pickBIC(object[as.character(G), modelNames, drop = FALSE], k = 3)
  if(all(is.na(bestBICs))) 
    { return(structure(NULL, bestBICvalues = bestBICs, prior = prior, 
                       control = control, initialization = initialization, 
                       class = "summary.mclustBIC"))
  }
  temp <- unlist(strsplit(names(bestBICs)[1], ","))
  bestModel <- temp[1] 
  G <- as.numeric(temp[2])
  if(G == 0) 
    { ans <- list(bic = bestBICs[1], 
                  z = unmap(rep(0,n)),
                  classification = rep(0, n), 
                  uncertainty = rep(0, n), 
                  n = n, d = ncol(data), 
                  modelName = bestModel, G = 0, 
                  loglik = n * logb(Vinv), 
                  Vinv = Vinv)
      orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "Vinv", 
                        "z", "classification", "uncertainty")
     return(structure(ans[orderedNames], bestBICvalues = bestBICs, 
                      prior = prior, control = control, 
                      initialization = initialization, 
                      class = "summary.mclustBIC"))
  }
  G1 <- G + 1
  z <- matrix(0, n, G1)
  if(d > 1 || !is.null(hcPairs))
    { z[!noise, 1:G] <- unmap(hclass(hcPairs, G)) }
  else 
    { z[!noise, 1:G] <- unmap(qclass(data[!noise], G)) }
  z[noise, G1] <- 1
  out <- me(modelName = bestModel, data = data, z = z, prior = prior, 
            control = control, warn = warn, Vinv = Vinv)
  obsNames <- if(is.null(dim(data))) 
                names(data) else dimnames(data)[[1]]
  classification <- map(out$z, warn = warn)
  classification[classification == G1] <- 0
  uncertainty <- 1 - apply(out$z, 1, max)
  names(classification) <- names(uncertainty) <- obsNames
  ans <- c(list(bic = as.vector(bestBICs[1]), classification = classification, 
                uncertainty = uncertainty, Vinv = Vinv), out)
  orderedNames <- c("modelName", "n", "d", "G", "bic", "loglik", "parameters", 
                    "Vinv", "z", "classification", "uncertainty")
  structure(ans[orderedNames], 
            bestBICvalues = bestBICs, 
            prior = prior, control = control, 
            initialization = initialization, 
            class = "summary.mclustBIC")
}

print.summary.mclustBIC <- function(x, digits = getOption("digits"), ...)
{
  if("classification" %in% names(x))
    { bic <- attr(x,"bestBICvalues")
      l <- length(bic)
      if(l == 1) 
        { cat("BIC value:\n")
          print(bic, digits = digits)
      }
      else {
          cat("Best BIC values:\n")
          bic <- drop(as.matrix(bic))
          bic <- rbind(BIC = bic, "BIC diff" = bic - max(bic))
          print(bic, digits = digits)
      }
      cat("\nClassification table for model (", colnames(bic)[1], "):", sep = "")
      print(table(x$classification), digits = digits, ...)
  }
  else {
    cat("Best BIC values:\n")
    x <- drop(as.matrix(x))
    x <- rbind(BIC = x, "BIC diff" = x - max(x))
    print(x, digits = digits)
  }
  invisible()
}

defaultPrior <- function(data, G, modelName, ...)
{
  aux <- list(...)
  if(is.null(aux$shrinkage)) {
    shrinkage <- 0.01
  }
  else if(is.na(aux$shrinkage) || !aux$shrinkage) {
    shrinkage <- 0
  }
  else if(aux$shrinkage < 0) {
    stop("negative value given for shrinkage")
  }
  else {
    shrinkage <- aux$shrinkage
  }
  if(is.null(aux$mean)) {
    mean <- if (is.null(dim(data))) 
      mean(data) else colMeans(data)
  }
  else if(any(is.na(aux$mean))) {
    if(shrinkage)
      stop("positive shrinkage with no prior mean specified")
    mean <- if (is.null(dim(data))) 
      mean(data) else colMeans(data)
  }
  else {
    if(!shrinkage)
      stop("prior mean specified but not shrinkage")
    mean <- aux$mean
  }
  switch(EXPR = modelName,
         E = ,
         V = ,
         X = {
           dof <- 3
           if(is.null(aux$scale)) {
             scale <- var(data)/G^2
           }
           else {
             scale <- aux$scale
           }
           list(shrinkage = shrinkage, mean = mean, dof = dof,
                scale = scale)
         },
         ##############################################################
         EII = ,
         VII = ,
         XII = ,
         EEI = ,
         EVI = ,
         VEI = ,
         VVI = ,
         XXI = {
           n <- nrow(data)
           p <- ncol(data)
           dof <- p + 2
           if(is.null(aux$scale)) {
             fac <- (1/G)^(2/p)
             scale <- (fac * sum(apply(data, 2, var)))/
               p
           }
           else {
             scale <- aux$scale
           }
           list(shrinkage = shrinkage, mean = mean, dof = dof,
                scale = scale)
         },
         ##############################################################
         EEE = ,
         EVE = ,
         VEE = ,
         VVE = ,
         EEV = ,
         VEV = ,
         EVV = ,
         VVV = ,
         XXX = {
           n <- nrow(data)
           p <- ncol(data)
           dof <- p + 2
           if(is.null(aux$scale)) {
             fac <- (1/G)^(2/p)
             if(n > p) {
               scale <- fac * var(data)
             }
             else {
               scale <- fac * diag(apply(data,
                                         2, var))
             }
           }
           else {
             scale <- aux$scale
           }
           list(shrinkage = shrinkage, mean = mean, dof = dof,
                scale = scale)
         },
         stop("no default prior for this model"))
}

emControl <- function(eps = .Machine$double.eps, 
                      tol = c(1.0e-05, sqrt(.Machine$double.eps)), 
                      itmax = c(.Machine$integer.max, .Machine$integer.max), 
                      equalPro = FALSE)
{
  if(any(eps < 0)) stop("eps is negative")
  if(any(eps >= 1))
    stop("eps is not less than 1")
  if(any(tol < 0))
    stop("tol is negative")
  if(any(tol >= 1))
    stop("tol is not less than 1")
  if(any(itmax < 0))
    stop("itmax is negative")
  if(length(tol) == 1)
    tol <- rep(tol, 2)
  if(length(itmax) == 1)
    itmax <- c(itmax, .Machine$integer.max)
  i <- is.infinite(itmax)
  if(any(i))
    itmax[i] <- .Machine$integer.max
  list(eps = eps, tol = tol, itmax = itmax, equalPro = equalPro)
}

priorControl <- function(functionName = "defaultPrior", ...)
{
  c(list(functionName = functionName), list(...))
}

cdensEEE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) 
  { WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "EEE", 
                     WARNING = WARNING, returnCode = 9))
  }
  if(is.null(parameters$variance$cholSigma))
    stop("variance parameters are missing")
  temp <- .Fortran("eseee",
                   as.logical(1),
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$cholSigma),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(p),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, G)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      if(warn) warning(WARNING)
    }
    z[] <- NA
    ret <- -9 
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EEE",
            WARNING = WARNING, retrunCode = ret)
}

emEEE <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepEEE(data, parameters = parameters, warn = warn)$z  
  meEEE(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepEEE <- function(data, parameters, warn = NULL, ...)
{
  if (is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) > 2)
    stop("data must be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EEE", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$cholSigma))
    stop("variance parameters are missing")
  temp <- .Fortran("eseee",
                   as.logical(1),
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$cholSigma),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(p),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      warning(WARNING)
      ret <- -4 
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      warning(WARNING)
      ret <- -5
    }
    z[] <- loglik <- NA
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "EEE", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meEEE <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEE", d = p, G = G, 
                     Sigma = matrix(as.double(NA), p, p), cholSigma = matrix(as.double(NA), p, p)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeee",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(p * p),
                     double(K),
                     double(p),
                     PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEE"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("meeeep",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(p * p),
                     double(K),
                     double(p),
                     PACKAGE = "mclust")[c(11:17, 10)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- matrix(temp[[6]], p, p)
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix( NA, p, p)
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- logprior <- NA
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(as.double(NA), p, p)
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    Sigma <- unchol(cholSigma, upper = TRUE)
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(Sigma) <- dimnames(cholSigma) <- 
    list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EEE", d = p, G = G,
                   sigma = sigma, Sigma = Sigma, cholSigma = cholSigma) 
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEE <- function(data, z, prior = NULL,  warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEE", d = p, G = G, 
                     sigma <- array(NA, c(p,p, G)), 
                     Sigma = matrix(as.double(NA), p, p), cholSigma = matrix(as.double(NA), p, p)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="EEE", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("mseee",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p),
                     double(p * G),
                     double(p * p),
                     double(G),
                     PACKAGE = "mclust")[7:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEE"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("mseeep",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     double(p),
                     double(p * G),
                     double(p * p),
                     double(G),
                     PACKAGE = "mclust")[11:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholSigma <- matrix(temp[[2]], p, p)
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  Sigma <- unchol(cholSigma, upper = TRUE)
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  WARNING <- NULL
  if(any(mu > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- sigma[] <- Sigma[] <- cholSigma[] <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(Sigma) <- dimnames(cholSigma) <- 
    list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EEE", d = p, G = G, 
                   sigma = sigma, Sigma = Sigma, cholSigma= cholSigma)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EEE", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
}

simEEE <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "EEE"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  if(is.null(cholSigma <- parameters$variance$cholSigma)) {
    if(is.null(Sigma <- parameters$variance$Sigma)) {
      stop("variance parameters must inlcude either Sigma or cholSigma"
      )
    }
    cholSigma <- chol(Sigma)
  }
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                              cholSigma, MARGIN = 2, STATS = mu[,k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "EEE")
}

cdensEEI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "EEI", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("eseei",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1 
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EEI",
            WARNING = WARNING, returnCode = ret)
}

cdensEII <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
  {
    if(is.null(warn)) warn <- mclust.options("warn")
    dimdat <- dim(data)
    if(is.null(dimdat) || length(dimdat) != 2)
      stop("data must be a matrix")
    data <- as.matrix(data)
    n <- nrow(data)
    p <- ncol(data)
    mu <- as.matrix(parameters$mean)
    G <- ncol(mu)
    if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
         any(is.null(parameters[c("pro", "mean", "variance")]))) {
      WARNING <- "parameters are missing"
      if(warn) warning(WARNING)
      z <- matrix(as.double(NA),n,G)
      dimnames(z) <- list(dimnames(data)[[1]], NULL)
      return(structure(z, logarithm = logarithm, modelName = "EII", 
                       WARNING = WARNING, returnCode = 9))
    }
    sigmasq <- parameters$variance$sigmasq
    if(sigmasq < 0)
      stop("sigma-squared is negative")
    if(!sigmasq) {
      WARNING <- "sigma-squared vanishes"
      if(warn) warning(WARNING)
      z <- matrix(as.double(NA),n,G)
      dimnames(z) <- list(dimnames(data)[[1]], NULL)
      return(structure(z, logarithm = logarithm, modelName = "EII", 
                       WARNING = WARNING, returnCode = 9))
    }
    temp <- .Fortran("eseii",
                     as.double(data),
                     as.double(mu),
                     as.double(sigmasq),
                     as.double(-1),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(-1),
                     double(1),
                     double(n * G),
                     PACKAGE = "mclust")[9:10]
    loglik <- temp[[1]]
    z <- matrix(temp[[2]], n, G)
    WARNING <- NULL
    if(loglik > signif(.Machine$double.xmax, 6)) {
      WARNING <- "sigma-squared falls below threshold"
      if(warn) warning(WARNING)
      z[] <- NA
      ret <- -1
    }
    else {
      if (!logarithm) z <- exp(z)
      ret <- 0
    }
    dimnames(z) <- list(dimnames(data)[[1]],NULL)
    structure(z, logarithm = logarithm, modelName = "EII",
              WARNING = WARNING, returnCode = ret)
  }

emEEI <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepEEI(data, parameters = parameters, warn = warn)$z  
  meEEI(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepEEI <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EEI", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("eseei",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "EEI", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meEEI <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data  should be in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEI", d = p, G = G, 
                     scale = NA, shape = rep(NA,p)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeei",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(p),
                     double(K),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEI"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("meeeip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(p),
                     double(K),
                     PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[6]]
  shape <- temp[[7]]
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(as.double(NA), p, p)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    Sigma <- matrix(as.double(NA), p, p)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EEI", d = p, G = G, 
                   sigma = sigma, Sigma = Sigma, 
                   scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEI", d = p, G = G, 
                     scale = NA, shape = rep(NA,p)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="EEI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
    
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("mseei",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p * G),
                     double(1),
                     double(p),
                     double(G),
                     PACKAGE = "mclust")[6:9]
  }
  else {
    storage.mode(z) <- "double"
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEI"), prior[names(
                                                          prior) != "functionName"]))
    temp <- .Fortran("mseeip",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(p * G),
                     double(1),
                     double(p),
                     double(G),
                     PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[2]]
  shape <- temp[[3]]
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(shape, scale) > signif(.Machine$double.xmax, 6)) || any(!c(
    scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- Sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else {
    sigma <- array(0, c(p, p, G))
    Sigma <- diag(scale * shape)
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EEI", d = p, G = G, 
                   sigma = sigma, Sigma = Sigma, 
                   scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EEI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
}

simEEI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "EEI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  shape <- parameters$variance$shape
  if(length(shape) != d)
    stop("shape incompatible with mean")
  cholSigma <- diag(sqrt(parameters$variance$scale * shape))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
                                       ncol = d) %*% cholSigma, MARGIN = 2, STATS = mu[, k],
                                FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "EEI")
}

cdensE <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
  mu <- drop(parameters$mean)
  G <- length(mu)
  if(any(is.na(unlist(parameters[c("mean", "variance")]))) ||
       any(is.null(parameters[c("mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(names(data), NULL)
    return(structure(z, logarithm = logarithm, modelName = "E",
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(length(sigmasq) > 1)
    if(warn) warning("more than one sigma-squared given")
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(names(data), NULL)
    return(structure(z, logarithm = logarithm, modelName = "E",
                     WARNING = WARNING, returnCode = 9))
  }
  temp <- .Fortran("es1e",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(-1),
                   as.integer(n),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  } 
  dimnames(z) <- list(names(data),NULL)
  structure(z, logarithm = logarithm, modelName = "E", 
            WARNING = WARNING, returnCode = ret) 
}

emE <- function(data, parameters, prior = NULL, control = emControl(), 
                warn = NULL, ...)
{
  z <- estepE(data, parameters = parameters, warn = warn)$z
  meE(data, z = z, prior = prior, control = control, 
      Vinv = parameters$Vinv, warn = warn)
}

estepE <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  
  mu <- drop(parameters$mean)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(length(sigmasq) > 1)
    if(warn) warning("more than one sigma-squared specified")
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "E", n=n, d=1, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("es1e",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(names(data),NULL) 
  structure(list(modelName = "E", n = n, d = 1, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

cdensEEV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "EEV", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape) ||
        is.null(parameters$variance$orientation)) 
    stop("variance parameters are missing")
  temp <- .Fortran("eseev",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(aperm(parameters$variance$orientation,c(2,1,3))),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(p),
                   double(p),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EEV",
            WARNING = WARNING, returnCode = ret)
}

emEEV <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepEEV(data, parameters = parameters, warn = warn)$z  
  meEEV(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepEEV <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EEV", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape) ||
        is.null(parameters$variance$orientation)) 
    stop("variance parameters are missing")
  temp <- .Fortran("eseev",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(aperm(parameters$variance$orientation,c(2,1,3))),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(p),
                   double(p),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1 
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "EEV", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meEEV <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEV", d = p, G = G, 
                     scale = NA, shape = rep(NA,p), orientation = array(NA,c(p,p,G)))
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
    
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p))
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeev",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     as.integer(lwork),
                     double(p * G),
                     double(1),
                     double(p),
                     double(p * p * G),
                     double(K),
                     double(lwork),
                     double(p),
                     PACKAGE = "mclust")[7:16]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEV"), 
                                                 prior[names(prior) !="functionName"]))
    temp <- .Fortran("meeevp",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     as.integer(lwork),
                     double(p * G),
                     double(1),
                     double(p),
                     double(p * p * G),
                     double(K),
                     double(lwork),
                     double(p),
                     PACKAGE = "mclust")[11:20]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  lapackSVDinfo <- temp[[5]]
  mu <- matrix(temp[[6]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[7]]
  shape <- temp[[8]]
  O <- aperm(array(temp[[9]], c(p, p, G)),c(2,1,3))
  pro <- temp[[10]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
    }
    z[] <- O[] <- shape[] <- NA
    scale <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    shape[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "a z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = FALSE)
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(O) <- list(dimnames(data)[[2]], dimnames(data)[[2]], 
                      NULL)
  ## Sigma = scale * O %*% diag(shape) %*% t(O)
  variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                   scale = scale, shape = shape, orientation = O) 
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv) 
  structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik),
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepEEV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EEV", d = p, G = G, 
                     scale = NA, shape = rep(NA,p), orientation=array(NA,c(p,p,G))) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="EEV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
    
  }
  #  shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
    "improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), G)
  if(is.null(prior)) {
    temp <- .Fortran("mseev",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(lwork),
                     as.integer(lwork),
                     double(p * G),
                     double(1),
                     double(p),
                     double(p * p * G),
                     double(G),
                     PACKAGE = "mclust")[7:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EEV"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("mseevp",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     double(lwork),
                     as.integer(lwork),
                     double(p * G),
                     double(1),
                     double(p),
                     double(p * p * G),
                     double(G),
                     PACKAGE = "mclust")[11:16]
  }
  lapackSVDinfo <- temp[[1]]
  mu <- matrix(temp[[2]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[3]]
  shape <- temp[[4]]
  O <- aperm( array(temp[[5]], c(p, p, G)), c(2,1,3))
  pro <- temp[[6]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
      if(warn) warning(WARNING)
      ret <- -4
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
      if(warn) warning(WARNING)
      ret <- -5
    }
    O[] <- shape[] <- scale <- NA
    sigma <- array(NA, c(p, p, G))
  }
  else if(any(c(abs(scale), shape) > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- scale <- O[] <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else {
    sigma <- scale * shapeO(shape, O, transpose = FALSE)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- dimnames(O) <- 
    list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
  variance <- list(modelName = "EEV", d = p, G = G, sigma = sigma,
                   scale = scale, shape = shape, orientation = O)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EEV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
}

simEEV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "EEV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  shape <- parameters$variance$shape
  if(length(shape) != d)
    stop("shape incompatible with mean")
  sss <- sqrt(parameters$variance$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    cholSigma <- t(parameters$variance$orientation[,  , k]) * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m,
                                       ncol = d) %*% cholSigma, MARGIN = 2, STATS = mu[, k],
                                FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "EEV")
}

emEII <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepEII(data, parameters = parameters, warn = warn)$z
  meEII(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepEII <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  p <- ncol(data)
  n <- nrow(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    if(warn) warning("variance parameters are missing")
  if(sigmasq < 0)
    stop("sigma-squared is negative")
  if(!sigmasq) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EII", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("eseii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "EII", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meEII <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  # number of groups
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meeii",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(K),
                     PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EII"), prior[names(prior) !=
                                                                                                 "functionName"]))
    temp <- .Fortran("meeiip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(K),
                     PACKAGE = "mclust")[c(11:17, 10)]
  }
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  Sigma <- diag(rep(sigmasq, p))
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
       sigmasq <= max(control$eps,0)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- Sigma
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL) 
  dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
  dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL) 
  variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                   Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance = variance, Vinv=Vinv)
  structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control, 
                 loglik = loglik),
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepEII <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EII", d = p, G = G, sigmasq = NA)
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="EII", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mseii",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p * G),
                     double(1),
                     double(G),
                     PACKAGE = "mclust")[6:8]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EII"), 
                                                 prior[names(prior) !="functionName"]))
    temp <- .Fortran("mseiip",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(p * G),
                     double(1),
                     double(G),
                     PACKAGE = "mclust")[10:12]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  Sigma <- diag(rep(sigmasq, p))
  for(k in 1:G)
    sigma[,  , k] <- Sigma
  WARNING <- NULL
  if(sigmasq > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]], NULL) 
  dimnames(mu) <- list(dimnames(data)[[2]], NULL) 
  dimnames(Sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]])
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL) 
  variance <- list(modelName = "EII", d = p, G = G, sigma = sigma, 
                   Sigma = Sigma, sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance = variance)
  structure(list(modelName = "EII", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
  
}

simEII <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d), modelName = "EII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  cholSigma <- diag(rep(sqrt(sigmasq), d))
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  cholSigma, MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "EII")
}

meE <- function(data, z, prior = NULL, control = emControl(), 
                Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be 1 dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "E", d = 1, G = G, sigmasq = NA)
    parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                          z=z, parameters=parameters, control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("me1e",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(G),
                     double(1),
                     double(K),
                     PACKAGE = "mclust")[6:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "E"), prior[names(prior) !=
                                                                                               "functionName"]))
    temp <- .Fortran("me1ep",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(G),
                     double(1),
                     double(K),
                     PACKAGE = "mclust")[c(10:16, 9)]
  }
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  ## log post <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
       sigmasq <= max(control$eps,0)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- logprior <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else if(its >= control$itmax[1]) {
    WARNING <- "iteration limit reached"
    if(warn) warning(WARNING)
    its <-  - its
    ret <- 1
  }
  else ret <- 0
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(names(data), NULL) 
  variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                 z = z, parameters = parameters, control = control, 
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
  
}

mstepE <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  # number of groups 
  G <- dimz[2]
  ##
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName="E", d=1, G=G, sigmasq=NA)
    parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                       variance=variance)
    return(structure(list(modelName="E", prior=prior, n=n, d=1, G=G,
                          z = z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("ms1e",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(G),
                     double(G),
                     double(1),
                     double(G),
                     PACKAGE = "mclust")[5:7]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "E"), prior[names(prior) !=
                                                                                               "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("ms1ep",
                     as.double(data),
                     z,
                     as.integer(n),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(G),
                     double(1),
                     double(G),
                     PACKAGE = "mclust")[9:11]
  }
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  WARNING <- NULL
  if(sigmasq > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    pro[] <- mu[] <- sigmasq <- NA
    ret <- -1
    
  }
  else ret <- 0
  dimnames(z) <- list(names(data), NULL) 
  variance <- list(modelName = "E", d = 1, G = G, sigmasq = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "E", prior = prior, n = n, d = 1, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
}

simE <- function(parameters, n, seed = NULL, ...)
{
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, 2), modelName = "E"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd)
  }
  structure(cbind(group = clabels, "1" = x), modelName = "E")
}

cdensEVI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- parameters$mean
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "EVI", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esevi",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "EVI",
            WARNING = WARNING, returnCode = ret)
}

emEVI <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepEVI(data, parameters = parameters, warn = warn)$z  
  meEVI(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepEVI <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "EVI", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esevi",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "EVI", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meEVI <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EVI", d = p, G = G, 
                     scale = NA, shape = matrix(as.double(NA),p,G)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("meevi",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EVI"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("meevip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(1),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      if(warn) warning("z column sum fell below threshold")
      WARNING <- "z column sum fell below threshold"
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EVI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control, 
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepEVI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "EVI", d = p, G = G, 
                     scale = NA, shape = matrix(as.double(NA),p,G)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="EVI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("msevi",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p * G),
                     double(1),
                     double(p * G),
                     double(G),
                     PACKAGE = "mclust")[6:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "EVI"), prior[names(
                                                          prior) != "functionName"]))
    temp <- .Fortran("msevip",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(p * G),
                     double(1),
                     double(p * G),
                     double(G),
                     PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!c(
    scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- scale <- shape[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1 
  }
  else {
    sigma <- array(apply(scale * shape, 2, diag), c(p, p, G))
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "EVI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "EVI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
}

simEVI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "EVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  shape <- as.matrix(parameters$variance$shape)
  if(!all(dim(shape) == dim(mean)))
    stop("shape incompatible with mean")
  sss <- sqrt(parameters$variance$scale * shape)
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  diag(sss[, k]), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "EVI")
}

mclust1Dplot <- function(data, parameters = NULL, z = NULL,
                         classification = NULL, truth = NULL, uncertainty = NULL, 
                         what = c("classification", "density", "errors", "uncertainty"), 
                         symbols = NULL, colors = NULL, ngrid = length(data),  
                         xlab = NULL,  xlim = NULL, CEX  = 1, 
                         main = FALSE, ...) 
{
  
  grid1 <- function (n, range = c(0, 1), edge = TRUE) 
  {
    if (any(n < 0 | round(n) != n)) 
      stop("n must be nonpositive and integer")
    G <- rep(0, n)
    if(edge) 
    { G <- seq(from = min(range), to = max(range), 
               by = abs(diff(range))/(n - 1)) }
    else 
    { lj <- abs(diff(range))
      incr <- lj/(2 * n)
      G <- seq(from = min(range) + incr, to = max(range) - incr, 
               by = 2 * incr) }
    G
  }
  
  densNuncer <- function(data, parameters) 
  {
    cden <- cdensV(data = data, parameters = parameters)
    if(parameters$variance$G != 1) 
    { z <- sweep(cden, MARGIN = 2, FUN = "*", STATS = parameters$pro)
      den <- apply(z, 1, sum)
      z <- sweep(z, MARGIN = 1, FUN = "/", STATS = den)
      data.frame(density = den, uncertainty = 1 - apply(z, 1, max))
    }
    else 
    { data.frame(density = cden, uncertainty =  rep(NA, length(cden))) }
  }
  
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if (is.null(xlab)) xlab <- " "
  p <- ncol(as.matrix(data))
  if (p != 1) 
    stop("for one-dimensional data only")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if (!is.null(parameters)) {
    mu <- parameters$mean
    L <- ncol(mu)
    sigmasq <- parameters$variance$sigmasq
    haveParams <- !is.null(mu) && !is.null(sigmasq) && 
      !any(is.na(mu)) && !any(is.na(sigmasq)) 
  }
  else haveParams <- FALSE
  if (is.null(xlim)) xlim <- range(data)
  if (haveParams) {
    G <- length(mu)
    if ((l <- length(sigmasq)) == 1) {
      sigmasq <- rep(sigmasq, G)
    }
    else if (l != G) {
      params <- FALSE
      warning("mu and sigma are incompatible")
    }
  }
  if (!is.null(truth)) {
    if (is.null(classification)) {
      classification <- truth
      truth <- NULL
    }
    else {
      if (length(unique(truth)) != 
            length(unique(classification))) truth <- NULL
      else truth <- as.character(truth)
    }
  }
  if(!is.null(classification)) 
  {
    classification <- as.character(classification)
    U <- sort(unique(classification))
    L <- length(U)
    if(is.null(symbols)) 
      { symbols <- rep("|", L) }
    else if(length(symbols) == 1) 
      { symbols <- rep(symbols, L) }
    else if(length(symbols) < L) 
      { warning("more symbols needed to show classification")
        symbols <- rep("|", L) }
    if(is.null(colors))
      { colors <- mclust.options("classPlotColors")[1:L] }
    else if(length(colors) == 1) 
      { colors <- rep(colors, L) }
    else if(length(colors) < L)
      { warning("more colors needed to show classification")
      colors <- rep("black", L) }
  }
  if (length(what) > 1) what <- what[1]
  choices <- c("classification", "density", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if (m) { 
    type <- choices[m] 
    bad <- what == "classification" && is.null(classification)
    bad <- bad || (what == "uncertainty" && is.null(uncertainty))
    bad <- bad || (what == "errors" && 
                     (is.null(classification) || is.null(truth)))
    if (bad) warning("insufficient input for specified plot")
  }
  else {
    bad <- !m
    warning("what improperly specified")
  }
  if (bad) what <- "bad"
  M <- L
  switch(EXPR = what,
         "classification" = 
         { plot(data, seq(from = 0, to = M, length = n), type = "n", 
                xlab = xlab, ylab = "", xlim = xlim, 
                ylim = grDevices::extendrange(r = c(0,M), f = 0.1),
                yaxt = "n", main = "", ...)
           axis(side = 2, at = 0:M, labels = c("", sort(unique(classification))))
           if(main) title("Classification")
           for(k in 1:L) 
              { I <- classification == U[k]
           points(data[I], rep(0, length(data[I])), 
                  pch = symbols[k], cex = CEX)
           points(data[I], rep(k, length(data[I])), 
                  pch = symbols[k], col = colors[k], cex = CEX)
           }
         },
         "errors" = 
         { ERRORS <- classError(classification, truth)$misclassified
           plot(data, seq(from = 0, to = M, length = n), type = "n", 
                xlab = xlab, ylab = "", xlim = xlim, yaxt = "n", main = "", ...)
           axis(side = 2, at = 0:M, labels = c("", unique(classification)))
           if(main) title("Classification Errors")
           good <- rep(TRUE, length(classification))
           good[ERRORS] <- FALSE
           sym <- "|"
           for(k in 1:L) 
           { K <- classification == U[k]
           I <- K & good
           if(any(I)) 
           { if(FALSE) 
           { sym <- if (L > 4) 
             1
           else if (k == 4) 
             5
           else k - 1
           }
             l <- sum(as.numeric(I))
             points(data[I], rep(0, l), pch = sym, 
                    col = colors[k], cex = CEX)
           }
           I <- K & !good
           if(any(I)) 
           { if(FALSE) 
           { sym <- if (L > 5) 
             16
           else k + 14 }
             l <- sum(as.numeric(I))
             points(data[I], rep(k, l), pch = sym, 
                    col = colors[k], cex = CEX)
             # points(data[I], rep(0, l), pch = sym, cex = CEX)
             # points(data[I], rep(-0.5, l), pch = sym, cex = CEX)
           }
           }
         },
         "uncertainty" = 
         { # x <- grid1(n = ngrid, range = xlim, edge = TRUE)
           # lx <- length(x)
           # Z <- densNuncer(data = x, parameters = parameters)
           # plot(x, Z$uncertainty, xlab = xlab, ylab = "uncertainty", 
           #      xlim = xlim, ylim = c(0,1), type = "l", main = "", ...)
           u <- (uncertainty - min(uncertainty))/
                (max(uncertainty) - min(uncertainty))
           b <- bubble(u, cex = CEX*c(0.3, 2), alpha = c(0.3, 1))
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data, uncertainty, type = "h", 
                xlab = xlab, ylab = "Uncertainty",
                xlim = xlim, ylim = c(0,1), main = "", 
                col = mapply(adjustcolor, 
                             col = colors[cl], 
                             alpha.f = b$alpha),
                ...)
           rug(data, lwd = 1, col = adjustcolor(par("fg"), alpha.f = 0.8))
           if(main) title("Uncertainty")
         },
         "density" = 
         { if(is.null(parameters$pro) && parameters$variance$G != 1) 
             stop("mixing proportions missing")
           x <- grid1(n = ngrid, range = xlim, edge = TRUE)
           lx <- length(x)
           Z <- densNuncer(data = x, parameters = parameters)
           plot(x, Z$density, xlab = xlab, ylab = "density", xlim = xlim, 
                type = "l", main = "", ...)
           if(main) title("Density")
         },
         { plot(data, rep(0, n), type = "n", xlab = "", ylab = "", 
                xlim = xlim, main = "", ...)
           points(data, rep(0, n), pch = "|", cex = CEX)
           if(main) title("Point Plot")
           # return(invisible())
         }
  )
  invisible()
}

mclust2Dplot <- function(data, parameters = NULL, z = NULL,
                         classification = NULL, truth = NULL, 
                         uncertainty = NULL, 
                         what = c("classification", "uncertainty", "errors"), 
                         addEllipses = TRUE,
                         symbols = NULL, colors = NULL, 
                         xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, 
                         scale = FALSE, CEX = 1, PCH = ".", 
                         main = FALSE, swapAxes = FALSE, ...) 
{
  if(dim(data)[2] != 2)
    stop("data must be two dimensional")
  if(is.null(classification) && !is.null(z))
    classification <- map(z)
  if(is.null(uncertainty) && !is.null(z))
    uncertainty <- 1 - apply(z, 1, max)
  if(!is.null(parameters)) 
    { mu <- parameters$mean
      L <- ncol(mu)
      sigma <- parameters$variance$sigma
      haveParams <- !is.null(mu) && !is.null(sigma) && 
                    !any(is.na(mu)) && !any(is.na(sigma)) 
  }
  else haveParams <- FALSE
  
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  if(is.null(xlim))
    xlim <- range(data[, 1])
  if(is.null(ylim))
    ylim <- range(data[, 2])
  if(scale) 
    { par(pty = "s")
      d <- diff(xlim) - diff(ylim)
      if(d > 0) { ylim <- c(ylim[1] - d/2, ylim[2] + d/2.) }
      else      { xlim <- c(xlim[1] + d/2, xlim[2] - d/2) }
  }
  
  dnames <- dimnames(data)[[2]]
  if(is.null(xlab)) 
    { xlab <- if(is.null(dnames)) "" else dnames[1] }
  if(is.null(ylab)) 
    { ylab <- if(is.null(dnames)) "" else dnames[2] }
  
  if(haveParams) 
    { G <- ncol(mu)
      dimpar <- dim(sigma)
      if(length(dimpar) != 3) 
        { haveParams <- FALSE
          warning("covariance must be a 3D matrix")
      }
      if(G != dimpar[3])
        { haveParams <- FALSE
          warning("means and variance parameters are incompatible")
      }
      mu <- array(mu, c(2, G))
      sigma <- array(sigma, c(2, 2, G))
    }
  
  if(swapAxes)
    { if(haveParams) 
        { mu <- mu[2:1,]
          sigma <- sigma[2:1, 2:1,]
      }
    data <- data[, 2:1]
  }
  
  if(!is.null(truth)) 
    { if(is.null(classification)) 
        { classification <- truth
          truth <- NULL
      }
      else 
        { if(length(unique(truth)) != 
             length(unique(classification))) 
                truth <- NULL
           else truth <- as.character(truth)
      }
  }
  
  if(charmatch("classification", what, nomatch = 0) && 
       is.null(classification) && !is.null(z))
    { classification <- map(z) }
  
  if(!is.null(classification)) 
    { classification <- as.character(classification)
      U <- sort(unique(classification))
      L <- length(U)
      noise <- (U[1] == "0")
      if(is.null(symbols))
      { if(L <= length(mclust.options("classPlotSymbols"))) 
      { symbols <- mclust.options("classPlotSymbols")[1:L]
        if(noise)
        { symbols <- c(16,symbols)[1:L] }
      }
      else if(L <= 9)
      { symbols <- as.character(1:9) }
      else if(L <= 26) { symbols <- LETTERS }
      }
      if(is.null(colors)) 
      { if(L <= length(mclust.options("classPlotColors"))) 
      { colors <- mclust.options("classPlotColors")[1:L]
        if(noise) 
        { colors <- unique(c("black", colors))[1:L] }
      }
      }
      else if(length(colors) == 1) colors <- rep(colors, L)
      if(length(symbols) < L) 
      { warning("more symbols needed to show classification ")
        symbols <- rep(16,L)
      }
      if(length(colors) < L) 
      { warning("more colors needed to show classification ")
        colors <- rep("black",L)
      }
  }
  
  if(length(what) > 1) what <- what[1]
  choices <- c("classification", "errors", "uncertainty")
  m <- charmatch(what, choices, nomatch = 0)
  if(m) 
    { what <- choices[m] 
      bad <- what == "classification" && is.null(classification)
      bad <- bad || (what == "uncertainty" && is.null(uncertainty))
      bad <- bad || (what == "errors" && 
                       (is.null(classification) || is.null(truth)))
      if(bad) warning("insufficient input for specified plot")
    }
    else 
    { bad <- !m
      warning("what improperly specified")
  }
  if(bad) what <- "bad"
  
  switch(EXPR = what,
         "classification" = {
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
           if(main) title("Classification")
           for(k in 1:L) 
              { I <- classification == U[k]
                points(data[I, 1], data[I, 2], pch = symbols[k], 
                       col = colors[k], cex = if(U[k] == "0") CEX/4 else CEX)
           }
         },
         "errors" = {
           ERRORS <- classError(classification, truth)$misclassified
           plot(data[, 1], data[, 2], type = "n", xlab = xlab, 
                ylab = ylab, xlim = xlim, ylim = ylim, main = "", ...)
           if(main) title("Classification Errors")
           CLASSES <- unique(as.character(truth))
           symOpen <- c(2, 0, 1, 5)
           symFill <- c(17, 15, 16, 18)
           good <- rep(TRUE,length(classification))
           good[ERRORS] <- FALSE
           if(L > 4) 
           { points(data[good, 1], data[good, 2], pch = 1, 
                    col = colors, cex = CEX)
             points(data[!good, 1], data[!good, 2], pch = 16, 
                    cex = CEX)
           }
           else 
           { for(k in 1:L) 
           { K <- truth == CLASSES[k]
             points(data[K, 1], data[K, 2], pch = symOpen[k], 
                    col = colors[k], cex = CEX)
             if(any(I <- (K & !good))) 
             { points(data[I, 1], data[I, 2], 
                      pch = symFill[k], cex = CEX)
             }
           }
           }
         },
         "uncertainty" = { 
           u <- (uncertainty - min(uncertainty))/
                (max(uncertainty) - min(uncertainty))
           b <- bubble(u, cex = CEX*c(0.3, 2), alpha = c(0.3, 0.9))
           cl <- sapply(classification, function(cl) which(cl == U))
           plot(data[, 1], data[, 2], pch = 19, 
                xlab = xlab, ylab = ylab, 
                xlim = xlim, ylim = ylim, main = "", 
                cex = b$cex, 
                col = mapply(adjustcolor, 
                             col = colors[cl], 
                             alpha.f = b$alpha),
                ...)
           if(main) title("Uncertainty")
         },
         {  plot(data[, 1], data[, 2], type = "n", 
                 xlab = xlab, ylab = ylab, 
                 xlim = xlim, ylim = ylim, main = "", ...)
            if(main) title("Point Plot")
            points(data[, 1], data[, 2], pch = PCH, cex = CEX)
         }
  )
  if(haveParams && addEllipses) 
    { ## plot ellipsoids
      for(k in 1:G) 
        mvn2plot(mu = mu[,k], sigma = sigma[,,k], k = 15)
  }

  invisible()
}

mvn2plot <- function(mu, sigma, k = 15, alone = FALSE, 
                     col = rep("grey30",3), pch = 8, lty = c(1,2), lwd = c(1,1)) 
{
  p <- length(mu)
  if (p != 2) 
    stop("only two-dimensional case is available")
  if (any(unique(dim(sigma)) != p)) 
    stop("mu and sigma are incompatible")
  ev <- eigen(sigma, symmetric = TRUE)
  s <- sqrt(rev(sort(ev$values)))
  V <- t(ev$vectors[, rev(order(ev$values))])
  theta <- (0:k) * (pi/(2 * k))
  x <- s[1] * cos(theta)
  y <- s[2] * sin(theta)
  xy <- cbind(c(x, -x, -x, x), c(y, y, -y, -y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
  if(alone) 
    { xymin <- apply(xy, 2, FUN = "min")
      xymax <- apply(xy, 2, FUN = "max")
      r <- ceiling(max(xymax - xymin)/2)
      xymid <- (xymin + xymax)/2
      plot(xy[, 1], xy[, 2], type = "n", xlab = "x", ylab = "y", 
           xlim = c(-r, r) + xymid[1], ylim = c(-r, r) + xymid[2])
  }
  l <- length(x)
  i <- 1:l
  for(k in 1:4) 
     { lines(xy[i,], col = col[1], lty = lty[1], lwd = lwd[1])
       i <- i + l
  }
  x <- s[1]
  y <- s[2]
  xy <- cbind(c(x, -x, 0, 0), c(0, 0, y, -y))
  xy <- xy %*% V
  xy <- sweep(xy, MARGIN = 2, STATS = mu, FUN = "+")
  lines(xy[1:2,], col = col[2], lty = lty[2], lwd = lwd[2])
  lines(xy[3:4,], col = col[2], lty = lty[2], lwd = lwd[2])
  points(mu[1], mu[2], col = col[3], pch = pch)
  invisible()
}

# new version
plot.Mclust <- function(x, 
                        what = c("BIC", "classification", "uncertainty", "density"), 
                        dimens = NULL, xlab = NULL, ylab = NULL, ylim = NULL,  
                        addEllipses = TRUE, main = TRUE,
                        ...) 
{
  
  object <- x # Argh.  Really want to use object anyway
  if(!inherits(object, "Mclust")) 
    stop("object not of class \"Mclust\"")
  
  data <- object$data
  p <- ncol(data)
  if(p == 1) 
    colnames(data) <- deparse(x$call$data)
  if(is.null(dimens)) 
    dimens <- seq(p)
  else
    dimens <- dimens[dimens <= p]
  d <- length(dimens)
  main <- if(is.null(main) || is.character(main)) FALSE else as.logical(main)
  
  ####################################################################
  what <- match.arg(what, several.ok = TRUE)
  oldpar <- par(no.readonly = TRUE)
  # on.exit(par(oldpar))
  
  plot.Mclust.bic <- function(...)
    plot.mclustBIC(object$BIC, xlab = xlab, ylim = ylim, ...)

  plot.Mclust.classification <- function(...)
  {  
    if(p == 1)
      { mclust1Dplot(data = data, 
                     # parameters = object$parameters, 
                     what = "classification",
                     classification = object$classification,
                     z = object$z, 
                     xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                     main = main, ...) 
    }
    if(p == 2) 
      { mclust2Dplot(data = data, what = "classification", 
                     classification = object$classification, 
                     parameters = if(addEllipses) object$parameters else NULL,
                     xlab = if(is.null(xlab)) colnames(data)[1] else xlab, 
                     ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                     main = main, ...) 
    }
    if(p > 2)
    { 
      if(d == 2)
      { if(addEllipses)
          { coordProj(data = data, what = "classification", 
                      parameters = object$parameters, z = object$z, 
                      dimens = dimens, main = main, ...) 
        }
        else
          { mclust2Dplot(data = data[,dimens], what = "classification", 
                         classification = object$classification,
                         # z = object$z, 
                         main = main, ...)
        }
      }
      else
      { if(addEllipses)
          { on.exit(par(oldpar))
            par(mfrow = c(d, d), 
                mar = rep(c(0.3,0.3/2),each=2), 
                oma = c(4, 4, 4, 4))
            for(i in seq(d))
            { for(j in seq(d)) 
              { if(i == j) 
                  { plot(data[,c(j,i)],type="n",xlab="",ylab="",axes=FALSE)
                    text(mean(par("usr")[1:2]),
                         mean(par("usr")[3:4]),
                         labels = colnames(data[,dimens])[i], 
                         cex=1.5, adj=0.5)
                    box()
                } 
                else 
                  { coordProj(data = data, 
                              what = "classification", 
                              parameters = object$parameters,
                              z = object$z,
                              dimens = dimens[c(j,i)], 
                              main = FALSE, 
                              xaxt = "n", yaxt = "n", ...)
                }
                if(i == 1 && (!(j%%2))) axis(3)
                if(i == d && (j%%2))    axis(1)
                if(j == 1 && (!(i%%2))) axis(2)
                if(j == d && (i%%2))    axis(4)
              }
            }
          }
      else
        { clPairs(data[,dimens], gap = 0.2, cex.labels = 1.5,
                  classification = object$classification, ...) }
      }
    }
    
  }

  plot.Mclust.uncertainty <- function(...) 
  {
    if(p == 1)
      { mclust1Dplot(data = data,
                     parameters = object$parameters,
                     z = object$z, what = "uncertainty", 
                     xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                     main = main, ...) 
    }
    if(p == 2) 
      { mclust2Dplot(data = data, parameters = object$parameters, 
                     z = object$z, what = "uncertainty", 
                     xlab = if(is.null(xlab)) colnames(data)[1] else xlab, 
                     ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                     addEllipses = addEllipses, main = main, ...)
    }
    if(p > 2) 
      { if(d == 2)
          { coordProj(data = data, parameters = object$parameters, 
                      z = object$z, what = "uncertainty", 
                      dimens = dimens, main = main, 
                      addEllipses = addEllipses, ...) }
        else
          { on.exit(par(oldpar))
            par(mfrow = c(d, d), 
                mar = rep(c(0.3,0.3/2),each=2), 
                oma = c(4, 4, 4, 4))
            for(i in seq(d))
               { for(j in seq(d)) 
                    { if(i == j) 
                        { plot(0,0,type="n",xlab="",ylab="",axes=FALSE)
                          text(0,0, colnames(data[,dimens])[i], 
                               cex=1.5, adj=0.5)
                          box()
                        } 
                      else 
                    { coordProj(data = data, 
                                what = "uncertainty", 
                                parameters = object$parameters, 
                                z = object$z,
                                dimens = dimens[c(j,i)], 
                                main = FALSE, 
                                addEllipses = addEllipses,
                                xaxt = "n", yaxt = "n", ...)
                    }
                    if(i == 1 && (!(j%%2))) axis(3)
                    if(i == d && (j%%2))   axis(1)
                    if(j == 1 && (!(i%%2))) axis(2)
                    if(j == d && (i%%2))   axis(4)
                  }
            }
          }
      }
  }
  
  plot.Mclust.density <- function(...)
  {
    if(p == 1)
      { mclust1Dplot(data = data,
                     parameters = object$parameters,
                     z = object$z, what = "density", 
                     xlab = if(is.null(xlab)) colnames(data)[dimens] else xlab, 
                     main = main, ...) 
    }
    if(p == 2) 
      { surfacePlot(data = data, parameters = object$parameters,
                    what = "density", nlevels = 11,
                    transformation = "log",
                    xlab = if(is.null(xlab)) colnames(data)[1] else xlab, 
                    ylab = if(is.null(ylab)) colnames(data)[2] else ylab,
                    main = main, ...) 
    }
    if(p > 2) 
      { objdens <- object
        objdens$varname <- colnames(data)
        objdens$range <- if(objdens$d > 1) apply(data, 2, range) else range(data)
        plotDensityMclustd(objdens, nlevels = 11, ...) }
  }
  
  if(interactive() & length(what) > 1)
    { title <- "Model-based clustering plots:"
      # present menu waiting user choice
      choice <- menu(what, graphics = FALSE, title = title)
      while(choice != 0)
           { if(what[choice] == "BIC")            plot.Mclust.bic(...)
             if(what[choice] == "classification") plot.Mclust.classification(...)
             if(what[choice] == "uncertainty")    plot.Mclust.uncertainty(...)
             if(what[choice] == "density")        plot.Mclust.density(...)
             # re-present menu waiting user choice
             choice <- menu(what, graphics = FALSE, title = title)
           }
  } 
  else 
    { if(any(what == "BIC"))            plot.Mclust.bic(...)
      if(any(what == "classification")) plot.Mclust.classification(...) 
      if(any(what == "uncertainty"))    plot.Mclust.uncertainty(...) 
      if(any(what == "density"))        plot.Mclust.density(...) 
  }
    
  invisible()
}

plot.mclustBIC <- function(x, G = NULL, modelNames = NULL, 
                           symbols = NULL, colors = NULL, 
                           xlab = NULL, ylab = "BIC", ylim = NULL, 
                           legendArgs = list(x = "bottomright", ncol = 2, cex = 1, inset = 0.01), 
                           ...)
{
  
  if(is.null(xlab)) xlab <- "Number of components"
  fill <- FALSE
  subset <- !is.null(attr(x, "initialization")$subset)
  noise <- !is.null(attr(x, "initialization")$noise)
  ret <- attr(x, "returnCodes") == -3
  legendArgsDefault <- list(x = "bottomright", ncol = 2, cex = 1, inset = 0.01)
  legendArgs <- append(as.list(legendArgs), legendArgsDefault)
  legendArgs <- legendArgs[!duplicated(names(legendArgs))]

  n <- ncol(x)
  dnx <- dimnames(x)
  x <- matrix(as.vector(x), ncol = n)
  dimnames(x) <- dnx
  if(is.null(modelNames))
    modelNames <- dimnames(x)[[2]]
  if(is.null(G))
    G <- as.numeric(dimnames(x)[[1]])
  # BIC <- x[as.character(G), modelNames, drop = FALSE]
  # X <- is.na(BIC)
  # nrowBIC <- nrow(BIC)
  # ncolBIC <- ncol(BIC)
  if(is.null(symbols)) 
  { colNames <- dimnames(x)[[2]]
    m <- length(modelNames)
    if(is.null(colNames)) 
    { symbols <- if(m > 9) LETTERS[1:m] else as.character(1:m)
      names(symbols) <- modelNames
    }
    else 
      { symbols <- mclust.options("bicPlotSymbols")[modelNames] }
  }
  if(is.null(colors)) 
  { colNames <- dimnames(x)[[2]]
    if(is.null(colNames)) 
    { colors <- 1:m
      names(colors) <- modelNames
    }
    else 
      { colors <- mclust.options("bicPlotColors")[modelNames] }
  }
  x <- x[,modelNames, drop = FALSE]
  if(is.null(ylim))
    ylim <- range(as.vector(x[!is.na(x)]))
  matplot(as.numeric(dnx[[1]]), x, type = "b", 
          xaxt = "n", xlim = range(G), ylim = ylim,
          pch = symbols, col = colors, lty = 1,
          xlab = xlab, ylab = ylab, main = "")
  axis(side = 1, at = as.numeric(dnx[[1]]))
  if(!is.null(legendArgs))
    { do.call("legend", c(list(legend = modelNames, col = colors, pch = symbols),
              legendArgs)) }
  invisible(symbols)
}


# old version: LS 20150317
sigma2decomp <- function(sigma, G = NULL, tol = sqrt(.Machine$double.eps), ...) 
{
  dimSigma <- dim(sigma)
  if(is.null(dimSigma)) 
    stop("sigma improperly specified")
  d <- dimSigma[1]
  if(dimSigma[2] != d) 
    stop("sigma improperly specified")
  l <- length(dimSigma)
  if(l < 2 || l > 3) 
    stop("sigma improperly specified")
  if(is.null(G)) 
    { if(l == 2) 
        { G <- 1
          sigma <- array(sigma, c(dimSigma, 1)) }
    else { G <- dimSigma[3] }
  } else 
    { if(l == 3 && G != dimSigma[3]) 
        stop("sigma and G are incompatible")
      if(l == 2 && G != 1) 
        sigma <- array(sigma, c(d,d,G))
  }
  
  # angle between subspaces
  subspace <- function(A, B)
  { for(k in 1:ncol(A))
       { B <- B - A[,k,drop=FALSE] %*% (t(A[,k,drop=FALSE]) %*% B) }
    norm(B, type = "2")
  }
  # check equality of values
  uniq <- function(x) { abs(max(x) - min(x)) < tol }
  
  decomp <- list(d = d, G = G, 
                 scale = rep(0, G), 
                 shape = matrix(0, d, G), 
                 orientation = array(0, c(d, d, G)))
  
  for(k in 1:G) 
     { ev <- eigen(sigma[,,k], symmetric = TRUE)
       temp <- log(ev$values); temp[!is.finite(temp)] <- 0
       logScale <- sum(temp)/d
       decomp$scale[k] <- exp(logScale)
       decomp$shape[,k] <- exp(temp - logScale)
       decomp$orientation[,,k] <- ev$vectors
  }
  scaleName <- "V"
  shapeName <- "V"
  orientName <- "V"
  # check scale/volume
  if(uniq(decomp$scale)) 
    { decomp$scale <- decomp$scale[1]
      scaleName <- "E"
  }
  # check shape
  if(all(apply(decomp$shape, 1, uniq))) 
    { decomp$shape <- decomp$shape[, 1]
      if(all(uniq(decomp$shape))) 
        { shapeName <- "I"
          decomp$shape <- rep(1, d)
      }
      else { shapeName <- "E" }
  }
  # check orientation
  eqOrientation <- 
  {  if(d == 2) all(apply(matrix(decomp$orientation, nrow = d * d, ncol = G), 
                                1, uniq))
    else       all(apply(decomp$orientation[,,-1,drop=FALSE], 3, 
                         function(o) subspace(decomp$orientation[,,1],o)) < tol)
  }
  if(eqOrientation)
    { decomp$orientation <- decomp$orientation[,,1]
      if(all(apply(cbind(decomp$orientation, diag(d)), 1, uniq)))
        { orientName <- "I"
          decomp$orientation <- NULL }
      else { orientName <- "E" }
  }
  
  decomp$modelName <- paste0(scaleName, shapeName, orientName)
  decomp$sigma <- sigma
  orderedNames <- c("sigma", "d", "modelName", "G", "scale", "shape", "orientation")
  return(decomp[orderedNames])
}

sigma2decomp <- function(sigma, G = NULL, tol = sqrt(.Machine$double.eps), ...) 
{
  dimSigma <- dim(sigma)
  if(is.null(dimSigma)) 
    stop("sigma improperly specified")
  d <- dimSigma[1]
  if(dimSigma[2] != d) 
    stop("sigma improperly specified")
  l <- length(dimSigma)
  if(l < 2 || l > 3) 
    stop("sigma improperly specified")
  if(is.null(G)) 
    { if(l == 2) 
        { G <- 1
          sigma <- array(sigma, c(dimSigma, 1)) }
    else { G <- dimSigma[3] }
  } else 
    { if(l == 3 && G != dimSigma[3]) 
        stop("sigma and G are incompatible")
      if(l == 2 && G != 1) 
        sigma <- array(sigma, c(d,d,G))
  }
  
  # angle between subspaces
  subspace <- function(A, B)
  { for(k in 1:ncol(A))
       { B <- B - A[,k,drop=FALSE] %*% (t(A[,k,drop=FALSE]) %*% B) }
    norm(B, type = "2")
  }
  # check equality of values
  uniq <- function(x) { abs(max(x) - min(x)) < tol }
  
  decomp <- list(d = d, G = G, 
                 scale = rep(0, G), 
                 shape = matrix(0, d, G), 
                 orientation = array(0, c(d, d, G)))
  
  for(k in 1:G) 
     { ev <- eigen(sigma[,,k], symmetric = TRUE)
       temp <- log(ev$values); temp[!is.finite(temp)] <- 0
       logScale <- sum(temp)/d
       decomp$scale[k] <- exp(logScale)
       decomp$shape[,k] <- exp(temp - logScale)
       decomp$orientation[,,k] <- ev$vectors
  }
  scaleName <- "V"
  shapeName <- "V"
  orientName <- "V"
  # check scale/volume
  if(uniq(decomp$scale)) 
    { decomp$scale <- decomp$scale[1]
      scaleName <- "E"
  }
  # check shape
  if(all(apply(decomp$shape, 1, uniq))) 
    { decomp$shape <- decomp$shape[, 1]
      if(all(uniq(decomp$shape))) 
        { shapeName <- "I"
          decomp$shape <- rep(1, d)
      }
      else { shapeName <- "E" }
  }
  # check orientation
  D <- decomp$orientation
  eqOrientation <- all(apply(D, 3, function(d) 
                       any(apply(d, 2, 
                                 function(x) cor(D[,,1], x)^2) > (1-tol))))
  if(eqOrientation)
    { decomp$orientation <- decomp$orientation[,,1]
      orientName <- "E"
      if(sum(abs(svd(decomp$orientation)$v) - diag(d)) < tol)
        { orientName <- "I"
          # decomp$orientation <- NULL 
        }
  }
  
  decomp$modelName <- paste0(scaleName, shapeName, orientName)
  decomp$sigma <- sigma
  orderedNames <- c("sigma", "d", "modelName", "G", "scale", "shape", "orientation")
  return(decomp[orderedNames])
}

decomp2sigma <- function(d, G, scale, shape, orientation = NULL, ...)
{
  nod <- missing(d)
  noG <- missing(G)
  lenScale <- length(scale)
  if(lenScale != 1) {
    if(!noG && G != lenScale)
      stop("scale incompatibile with G")
    G <- lenScale
  }
  shape <- as.matrix(shape)
  p <- nrow(shape)
  if(!nod && p != d)
    stop("shape incompatible with d")
  d <- p
  g <- ncol(shape)
  if(g != 1) {
    if(!is.null(G) && g != G)
      stop("shape incompatible with scale")
    if(!noG && g != G)
      stop("shape incompatible with G")
    G <- g
  }
  if(is.null(orientation)) {
    orientName <- "I"
    if(is.null(G)) {
      G <- if(noG) 1 else G
    }
    orientation <- array(diag(d), c(d, d, G))
  }
  else {
    dimO <- dim(orientation)
    l <- length(dimO)
    if(is.null(dimO) || l < 2 || l > 3 || dimO[1] != dimO[2])
      stop("orientation improperly specified")
    if(dimO[1] != d)
      stop("orientation incompatible with shape")
    if(l == 3) {
      orientName <- "V"
      if(is.null(G)) {
        if(!noG && dimO[3] != G)
          stop("orientation incompatible with G")
        G <- dimO[3]
      }
      else if(G != dimO[3])
        stop("orientation incompatible with scale and/or shape"
        )
    }
    else {
      orientName <- "E"
      if(is.null(G)) {
        G <- if(noG) 1 else G
      }
      orientation <- array(orientation, c(d, d, G))
    }
  }
  if(G == 1) {
    scaleName <- shapeName <- "X"
  }
  else {
    scaleName <- if(lenScale == 1) "E" else "V"
    shapeName <- if(g == 1) "E" else "V"
    scale <- rep(scale, G)
    shape <- matrix(shape, nrow = d, ncol = G)
  }
  sigma <- array(0, c(d, d, G))
  for(k in 1:G) 
     { sigma[,,k] <- crossprod(t(orientation[,,k]) * sqrt(scale[k] * shape[,k])) }

  structure(sigma, modelName = paste0(scaleName, shapeName, orientName))
}

grid1 <- function (n, range = c(0, 1), edge = TRUE) 
{
  if (any(n < 0 | round(n) != n)) 
    stop("n must be nonpositive and integer")
  G <- rep(0, n)
  if (edge) {
    G <- seq(from = min(range), to = max(range), by = abs(diff(range))/(n - 
                                                                          1))
  }
  else {
    lj <- abs(diff(range))
    incr <- lj/(2 * n)
    G <- seq(from = min(range) + incr, to = max(range) - 
               incr, by = 2 * incr)
  }
  G
}

grid2 <- function (x, y) 
{
  lx <- length(x)
  ly <- length(y)
  xy <- matrix(0, nrow = lx * ly, ncol = 2)
  l <- 0
  for (j in 1:ly) {
    for (i in 1:lx) {
      l <- l + 1
      xy[l,] <- c(x[i], y[j])
    }
  }
  xy
}

hypvol <- function (data, reciprocal = FALSE) 
{
  dimdat <- dim(data)
  oneD <- (is.null(dimdat) || length(dimdat[dimdat > 1]) == 1)
  if (oneD) {
    n <- length(as.vector(data))
    if (reciprocal) {
      ans <- 1/diff(range(data))
    }
    else {
      ans <- diff(range(data))
    }
    return(ans)
  }
  if (length(dimdat) != 2) 
    stop("data must be a vector or a matrix")
  data <- as.matrix(data)
  
  sumlogdifcol <- function(x) 
    sum(log(apply(x, 2, function(colm) diff(range(colm)))))
  
  bdvolog <- sumlogdifcol(data)
  pcvolog <- sumlogdifcol(princomp(data)$scores)
  
  volog <- min(bdvolog,pcvolog)
  
  if (reciprocal) {
    minlog <- log(.Machine$double.xmin)
    if (-volog < minlog) {
      warning("hypervolume smaller than smallest machine representable positive number")
      ans <- 0
    }
    else ans <- exp(-volog)
  }
  else {
    maxlog <- log(.Machine$double.xmax)
    if (volog > maxlog) {
      warning("hypervolume greater than largest machine representable number")
      ans <- Inf
    }
    else ans <- exp(volog)
  }
  
  ans
}

imputeData <- function(x, categorical=NULL, seed=NULL) 
{
  if(!requireNamespace("mix", quietly = TRUE)) 
     stop("imputeData function require 'mix' package to be installed!")
  
  fac <- apply( x, 2, is.factor)
  if(is.null(categorical)) 
    { categorical <- fac }
  else 
    { if(any(!categorical & fac)) 
        { stop("x has a factor that is not designated as categorical") }
      if(any(categorical | !fac)) 
        { warning("a categorical is not designated as a factor")
         for(i in which(categorical | !fac)) 
             x[[i]] <- as.factor(x[[i]])
      }
  }
  
  # remove categorical variables and add dummy variable
  if(nocat <- !any(categorical)) 
    { x <- cbind(as.factor(1),x)
      categorical <- c(TRUE, categorical)
  }
  
  ord <- c(which(categorical),which(!categorical))
  
  # do the imputations
  s <- mix::prelim.mix(x[,ord],p=sum(categorical))
  if(is.null(seed)) seed <- runif(1, min = .Machine$integer.max/1024,
                                     max = .Machine$integer.max)
  mix::rngseed(seed) # set random number generator seed
  thetahat <- mix::em.mix(s) # find ML estimate
  newtheta <- mix::da.mix(s, thetahat, steps=100, showits=TRUE)
  ximp <- mix::imp.mix(s, newtheta) # impute under newtheta
  
  if(nocat) ximp[,-1] else ximp[,order(ord)]
}

"[.mclustBIC" <- function (x, i, j, drop = FALSE) 
{
  ATTR <- attributes(x)[c("G", "modelNames", "prior", "control", 
                          "initialization", "Vinv", "warn", "n", "d", 
                          "oneD", "returnCodes", "class")]
  oldClass(x) <- NULL
  x <- NextMethod("[")
  if (is.null(dim(x))) return(x)
  ATTR$G <- as.numeric(dimnames(x)[[1]])
  ATTR$modelNames <- dimnames(x)[[2]]
  ATTR$returnCodes <- ATTR$returnCodes[dimnames(x)[[1]],dimnames(x)[[2]],
                                       drop=FALSE]
  do.call("structure", c(list(.Data = x), ATTR))
}

bic <- function(modelName, loglik, n, d, G, noise = FALSE, equalPro = FALSE, ...)
{
  nparams <- nMclustParams(modelName = modelName, d = d, G = G, 
                           noise = noise, equalPro = equalPro)
  2 * loglik - nparams * log(n)
}

checkModelName <- function(modelName)
{
  switch(EXPR = modelName,
         "X" = ,
         "E" = ,
         "V" = ,
         "XII" = ,
         "XXI" = ,
         "XXX" = ,
         "EII" = ,
         "VII" = ,
         "EEI" = ,
         "VEI" = ,
         "EVI" = ,
         "VVI" = ,
         "EEE" = ,
         "EVE" = ,
         "VEE" = ,
         "VVE" = ,
         "EEV" = ,
         "VEV" = ,
         "EVV" = ,
         "VVV" = TRUE,
         stop("invalid model name"))
}

em <- function(modelName, data, parameters, prior = NULL, control = emControl(), 
               warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("em", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

estep <- function(modelName, data, parameters, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("estep", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}


#############################################################################

mclustVariance <- function(modelName, d=NULL, G=2)
{
  x <- -1
  if (nchar(modelName) == 1) {
    if (!is.null(d) && d != 1)  stop("modelName and d are incompatible")
    varList <- switch(EXPR = modelName,
                      "X" = list(sigmasq = x),
                      "E" = list(sigmasq = x),
                      "V" = list(sigmasq = rep(x,G)),
                      stop("modelName not recognized"))
  }
  else {
    if (nchar(modelName) != 3) stop("modelName is misspecified")
    if (is.null(d)) d <- 3
    varList <- switch(EXPR = modelName,
                      "XII" = list(sigmasq = x),
                      "EII" = list(sigmasq = x, 
                                   scale = x, 
                                   shape = rep(x,d)),
                      "VII" = list(sigmasq = rep(x,G),
                                   scale = rep(x,G),
                                   shape = rep(x,d)),
                      "XXI" = list(scale = x, 
                                   shape = rep(x,d)),
                      "EEI" = list(scale = x, 
                                   shape = rep(x,d)),
                      "EVI" = list(scale = x, 
                                   shape = matrix(x,d,G)),
                      "VEI" = list(scale = rep(x,G), 
                                   shape = rep(x,d)),
                      "VVI" = list(scale = rep(x,G), 
                                   shape = matrix(x,d,G)),
                      "XXX" = { M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                                list(cholSigma = M) },
                      "EEE" = { M <- matrix(x,d,d); M[row(M) > col(M)] <- 0;
                                list(cholSigma = M) },
                      "VEE" = list(scale = rep(x,G), 
                                   shape = rep(x,d),
                                   orientation = matrix(x,d,d)),
                      "VVE" = list(scale = rep(x,G), 
                                   shape = matrix(x,d,G),
                                   orientation = matrix(x,d,d)),
                      "EVV" = list(scale = x, 
                                   shape = matrix(x,d,G),
                                   orientation = array(x,c(d,d,G))),
                      "EVE" = list(scale = x, 
                                   shape = matrix(x,d,G),
                                   orientation = matrix(x,d,d)),
                      "EEV" = list(scale = x, 
                                   shape = rep(x,d),
                                   orientation = array(x,c(d,d,G))),
                      "VEV" = list(scale = x, 
                                   shape = matrix(x,d,G),
                                   orientation = array(x,c(d,d,G))),
                      "VVV" = { A <- array(x,c(d,d,G));
                                I <- row(A[,,1]) > col(A[,,1])
                                for (k in 1:G) A[,,k][I] <- 0
                                list(cholsigma = A)},
                      stop("modelName not recognized"))
  }
  c(modelName = modelName, d = d, G = G, varList)
}

me <- function(modelName, data, z, prior = NULL, control = emControl(), 
               Vinv = NULL, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("me", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

mstep <- function(modelName, data, z, prior = NULL, warn = NULL, ...)
{
  checkModelName(modelName)
  funcName <- paste("mstep", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

mvn <- function(modelName, data, prior = NULL, warn = NULL, ...)
{
  modelName <- switch(EXPR = modelName,
                      "E" = "X",
                      "V" = "X",
                      "X"  =  "X",
                      "Spherical" = "XII",
                      "EII" = "XII",
                      "VII" = "XII",
                      "XII" = "XII",
                      "Diagonal" = "XXI",
                      "EEI" = "XXI",
                      "VEI" = "XXI",
                      "EVI" = "XXI",
                      "VVI" = "XXI",
                      "XXI" = "XXI",
                      "Ellipsoidal" = "XXX",
                      "EEE" = "XXX",
                      "VEE" = "XXX",
                      "EVE" = "XXX",
                      "EVV" = "XXX",
                      "VVE" = "XXX",
                      "EEV" = "XXX",
                      "VEV" = "XXX",
                      "VVV" = "XXX",
                      "XXX" = "XXX",
                      stop("invalid model name"))

  funcName <- paste("mvn", modelName, sep = "")
  mc <- match.call()
  mc[[1]] <- as.name(funcName)
  mc[[2]] <- NULL
  out <- eval(mc, parent.frame())
  varnames <- colnames(as.matrix(data))
  if(!all(is.null(varnames)))
  { rownames(out$parameters$mean) <- varnames
    dimnames(out$parameters$variance$Sigma) <- list(varnames, varnames)
    dimnames(out$parameters$variance$sigma) <- list(varnames, varnames, NULL)
  }
  return(out)
}

nVarParams <- function(modelName, d, G, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  # checkModelName(modelName)
  switch(EXPR = modelName,
         "E"   = 1,
         "V"   = G,
         "EII" = 1,
         "VII" = G,
         "EEI" = d,
         "VEI" = G + (d-1),
         "EVI" = 1 + G * (d-1),
         "VVI" = G * d,
         "EEE" = d*(d+1)/2,
         "EVE" = 1 + G*(d-1) + d*(d-1)/2,
         "VEE" = G + (d-1) + d*(d-1)/2,
         "VVE" = G + G * (d-1) + d*(d-1)/2,
         "EEV" = 1 + (d-1) + G * d*(d-1)/2,
         "VEV" = G + (d-1) + G * d*(d-1)/2,
         "EVV" = 1 - G + G * d*(d+1)/2,
         "VVV" = G * d*(d+1)/2,
         stop("invalid model name"))
}

nMclustParams <- function(modelName, d, G, noise = FALSE, equalPro = FALSE, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  if(G == 0) 
    { ## one noise cluster case
      if(!noise) stop("undefined model")
      nparams <- 1
  }
  else 
    { nparams <- nVarParams(modelName, d = d, G = G) + G*d
      if(!equalPro)
        nparams <- nparams + (G - 1)
      if(noise)
        nparams <- nparams + 2
  }
  return(nparams)
}

sim <- function(modelName, parameters, n, seed = NULL, ...)
{
  modelName <- switch(EXPR = modelName,
                      X = "E",
                      XII = "EII",
                      XXI = "EEI",
                      XXX = "EEE",
                      modelName)
  checkModelName(modelName)
  funcName <- paste("sim", modelName, sep = "")
  mc <- match.call(expand.dots = TRUE)
  mc[[1]] <- as.name(funcName)
  mc$modelName <- NULL
  eval(mc, parent.frame())
}

cdensVEI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "VEI", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvei",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],dimnames(mu)[[2]])
  structure(z, logarithm = logarithm, modelName = "VEI",
            WARNING = WARNING, returnCode = ret)
}

emVEI <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepVEI(data, parameters = parameters, warn = warn)$z  
  meVEI(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepVEI <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  } 
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VEI", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvei",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "VEI", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meVEI <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VEI", d = p, G = G, 
                     scale = rep(NA,G), shape = rep(NA,p)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevei",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax),
                     as.double(control$tol),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p),
                     double(K),
                     double(G),
                     double(p),
                     double(p * G),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = data, 
                                                      G = G, 
                                                      modelName = "VEI"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("meveip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax),
                     as.double(control$tol),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p),
                     double(K),
                     double(G),
                     double(p),
                     double(p * G),
                     PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]][1]
  inner <- temp[[2]][2]
  err <- temp[[3]][1]
  inerr <- temp[[3]][2]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- temp[[7]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= control$itmax[2]) {
      WARNING <- "inner iteration limit reached"
      if(warn) warning(WARNING)
      inner <-  - inner
      ret <- 2
    }
    else if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  attr(info, "inner") <- c(iterations = inner, error = inerr)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
  variance <- list(modelName = "VEI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepVEI <- function(data, z, prior = NULL, warn = NULL, control = NULL,...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VEI", d = p, G = G, 
                     scale = rep(NA,G), shape = rep(NA,p)) 
    parameters <- list(pro=rep(NA,G), 
                       mean=matrix(as.double(NA),p,G), variance=variance)
    return(structure(list(modelName="VEI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
    
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if (is.null(control)) control <- emControl()
  itmax <- if(length(control$itmax) == 1) control$itmax else control$
    itmax[2]
  tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
  if(is.null(prior)) {
    temp <- .Fortran("msvei",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.integer(itmax),
                     as.double(tol),
                     double(p * G),
                     double(G),
                     double(p),
                     double(G),
                     double(G),
                     double(p),
                     double(p * G),
                     PACKAGE = "mclust")[6:11]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VEI"), prior[names(
                                                          prior) != "functionName"]))
    temp <- .Fortran("msveip",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     as.integer(itmax),
                     as.double(tol),
                     double(p * G),
                     double(G),
                     double(p),
                     double(G),
                     double(G),
                     double(p),
                     double(p * G),
                     PACKAGE = "mclust")[10:15]
  }
  inner <- temp[[1]]
  inerr <- temp[[2]]
  mu <- matrix(temp[[3]], p, G)
  scale <- temp[[4]]
  shape <- temp[[5]]
  dimnames(mu) <- list(NULL, as.character(1:G))
  pro <- temp[[6]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
                                                                     c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else {
    ret <- 0
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(scale[k] * shape)
    if(inner >= itmax) {
      WARNING <- "inner iteration limit reached"
      if(warn) warning(WARNING)
      inner <-  - inner
    }
  }
  info <- c(iterations = inner, error = inerr)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VEI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VEI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control),
            info = info, WARNING = WARNING, returnCode = ret)
  
}

simVEI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "VEI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  diag(rtscale[k] * rtshape), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VEI")
}

cdensV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
  mu <- drop(parameters$mean)
  G <- length(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
     || any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(names(data), NULL)
    return(structure(z, logarithm = logarithm, modelName = "V", 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(names(data), NULL)
    return(structure(z, logarithm = logarithm, modelName = "V", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (length(sigmasq) == 1) sigmasq <- rep(sigmasq,G)
  temp <- .Fortran("es1v",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(-1),
                   as.integer(n),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(names(data),NULL)
  structure(z, logarithm = logarithm, modelName = "V",
            WARNING = WARNING, returnCode = ret)
}

emV <- function(data, parameters, prior = NULL, control = emControl(), 
                warn = NULL, ...)
{
  z <- estepV(data, parameters = parameters, warn = warn)$z  
  meV(data, z = z, prior = prior, control = control, 
      Vinv = parameters$Vinv, warn = warn)
}

estepV <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- drop(data)
  n <- length(data)
  pro <- parameters$pro 
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- drop(parameters$mean)
  G <- length(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv 
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")])))
     || any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(names(data), NULL)
    return(structure(list(modelName = "V", n=n, d=1, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("es1v",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[8:9]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) 
    { WARNING <- "cannot compute E-step"
      if(warn) warning(WARNING)
      z[] <- loglik <- NA
      ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(names(data),NULL)
  structure(list(modelName = "V", n = n, d = 1, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

cdensVEV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "VEV", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape) ||
        is.null(parameters$variance$orientation)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvev",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(aperm(parameters$variance$orientation,c(2,1,3))),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(p),
                   double(p),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL) 
  structure(z, logarithm = logarithm, modelName = "VEV",
            WARNING = WARNING, returnCode = ret)
}

emVEV <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepVEV(data, parameters = parameters, warn = warn)$z  
  meVEV(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepVEV <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VEV", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape) ||
        is.null(parameters$variance$orientation)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvev",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(aperm(parameters$variance$orientation,c(2,1,3))),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(p),
                   double(p),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[13:14]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1  
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL) 
  structure(list(modelName = "VEV", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
  
}

meVEV <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VEV", d = p, G = G, 
                     scale=rep(NA,G), shape=rep(NA,p), orientation=array(NA,c(p,p,G))) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters,
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevev",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax),
                     as.double(control$tol),
                     as.double(control$eps),
                     as.integer(lwork),
                     double(p * G),
                     double(G),
                     double(p),
                     double(p * p * G),
                     double(K),
                     double(lwork),
                     double(p),
                     PACKAGE = "mclust")[7:16]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VEV"), prior[names(prior) !=
                                                                                                 "functionName"]))
    temp <- .Fortran("mevevp",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax),
                     as.double(control$tol),
                     as.double(control$eps),
                     as.integer(lwork),
                     double(p * G),
                     double(G),
                     double(p),
                     double(p * p * G),
                     double(K),
                     double(lwork),
                     double(p),
                     PACKAGE = "mclust")[11:20]
  }
  z <- temp[[1]]
  its <- temp[[2]][1]
  inner <- temp[[2]][2]
  err <- temp[[3]][1]
  inerr <- temp[[3]][2]
  loglik <- temp[[4]]
  lapackSVDinfo <- temp[[5]]
  mu <- matrix(temp[[6]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[7]]
  shape <- temp[[8]]
  O <- aperm( array(temp[[9]], c(p, p, G)), c(2,1,3))
  pro <- temp[[10]]
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
    }
    if(warn) warning(WARNING)
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    O[] <- shape[] <- scale[] <- NA
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- shapeO(shape, O, transpose = FALSE)
    sigma <- sweep(sigma, MARGIN = 3, STATS = scale, FUN = "*")
    if(inner >= control$itmax[2]) {
      WARNING <- "inner iteration limit reached"
      if(warn) warning(WARNING)
      inner <-  - inner
      ret <- 2
    }
    else if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- structure(c(iterations = its, error = err), inner = c(
    iterations = inner, error = inerr))
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- dimnames(O) <- 
    list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
  ##  Sigma = scale * O %*% diag(shape) %*% t(O)
  variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                   scale = scale, shape = shape, orientation = O)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv) 
  structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepVEV <- function(data, z, prior = NULL, warn = NULL, control = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VEV", d = p, G = G, 
                     scale = rep(NA,G), shape = rep(NA,p), orientation = array(NA,c(p,p,G))) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="VEV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
    
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    return(structure(list(n = n, d = p, G = G, mu = matrix(as.double(NA),
                                                           p, G), sigma = array(NA, c(p, p, G)), decomp = list(
                                                             d = p, G = G, scale = rep(NA, G), shape = rep(NA, p),
                                                             orientation = array(NA, c(p, p, G))), pro = rep(NA,
                                                                                                             G), modelName = "VEV", prior = prior), WARNING = 
                       WARNING))
  }
  #  shape <- sqrt(rev(sort(shape/exp(sum(log(shape))/p))))
  if(any(is.na(z)) || any(z < 0) || any(z > 1)) stop(
    "improper specification of z")
  if (is.null(control)) control <- emControl()
  itmax <- if(length(control$itmax) == 1) control$itmax else control$
    itmax[2]
  tol <- if(length(control$tol) == 1) control$tol else control$tol[2]
  lwork <- max(3 * min(n, p) + max(n, p), 5 * min(n, p), p + G)
  if(is.null(prior)) {
    temp <- .Fortran("msvev",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(lwork),
                     as.integer(lwork),
                     as.integer(itmax),
                     as.double(tol),
                     double(p * G),
                     double(G),
                     double(p),
                     double(p * p * G),
                     double(G),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VEV"), prior[names(prior) !=
                                                                                                 "functionName"]))
    temp <- .Fortran("msvevp",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     double(lwork),
                     as.integer(lwork),
                     as.integer(itmax),
                     as.double(tol),
                     double(p * G),
                     double(G),
                     double(p),
                     double(p * p * G),
                     double(G),
                     PACKAGE = "mclust")[11:18]
  }
  lapackSVDinfo <- temp[[1]]
  inner <- temp[[2]]
  inerr <- temp[[3]]
  mu <- matrix(temp[[4]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  scale <- temp[[5]]
  shape <- temp[[6]]
  O <- aperm(array(temp[[7]], c(p, p, G)),c(2,1,3))
  pro <- temp[[8]] 
  WARNING <- NULL
  if(lapackSVDinfo) {
    if(lapackSVDinfo > 0) {
      WARNING <- "LAPACK DGESVD fails to converge"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DGESVD"
      if(warn) warning(WARNING)
    }
    O[] <- shape[] <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -9
  }
  else if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(
    !c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- O[] <- shape[] <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else {
    sigma <- sweep(shapeO(shape, O, transpose = FALSE), MARGIN = 3,
                   STATS = scale, FUN = "*")
    if(inner >= itmax) {
      WARNING <- "inner iteration limit reached"
      if(warn) warning(WARNING)
      inner <-  - inner
    }
    ret <- 2
  }
  info <- c(iteration = inner, error = inerr)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- dimnames(O) <-
    list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
  variance <- list(modelName = "VEV", d = p, G = G, sigma = sigma, 
                   scale = scale, shape = shape, orientation = O)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VEV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control),
            info = info, WARNING = WARNING, returnCode = ret)
  
}

simVEV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "VEV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(length(rtshape) != d)
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    sss <- rtscale[k] * rtshape
    cholSigma <- t(parameters$variance$orientation[,  , k]) * sss
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  cholSigma, MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VEV")
}

cdensVII <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "VII", 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq 
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "VII", 
                     WARNING = WARNING, returnCode = 9))
  }
  temp <- .Fortran("esvii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0 
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "VII",
            WARNING = WARNING, returnCode = ret)
}

emVII <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepVII(data, parameters = parameters, warn = warn)$z  
  meVII(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepVII <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  sigmasq <- parameters$variance$sigmasq 
  if(is.null(sigmasq))
    stop("variance parameters are missing")
  if(any(sigmasq < 0))
    stop("sigma-squared is negative")
  if(any(!sigmasq)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VII", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = -1))
  }
  temp <- .Fortran("esvii",
                   as.double(data),
                   as.double(mu),
                   as.double(sigmasq),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[9:10]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "VII", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meVII <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data must be in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevii",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(K),
                     PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VII"), prior[names(prior) !=
                                                                                                 "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("meviip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(K),
                     PACKAGE = "mclust")[c(11:17, 10)]
  }
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
       any(sigmasq <= max(control$eps, 0))) {
    WARNING <- "sigma-squared falls below threshold"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(0, c(p, p, G))
    for(k in 1:G)
      sigma[,  , k] <- diag(rep(sigmasq[k], p))
    if(its >= control$itmax[1]) 
      { WARNING <- "iteration limit reached"
        if(warn) warning(WARNING)
        its <-  - its
        ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VII", d = p, G = G, 
                   sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control, 
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

meVVI <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VVI", d = p, G = G, 
                     scale = rep(NA,G), shape = matrix(as.double(NA),p,G)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvi",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VVI"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VVI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepVII <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal number of observations")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("msvii",
                     as.double(data),
                     z,
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p * G),
                     double(G),
                     double(G),
                     PACKAGE = "mclust")[6:8]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VII"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("msviip",
                     as.double(data),
                     z,
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(p * G),
                     double(G),
                     double(G),
                     PACKAGE = "mclust")[10:12]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  sigma <- array(0, c(p, p, G))
  for(k in 1:G)
    sigma[,  , k] <- diag(rep(sigmasq[k], p))
  WARNING <- NULL
  if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VII", d = p, G = G, 
                   sigma = sigma, sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VII", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters),
            WARNING = WARNING, returnCode = ret)
}

simVII <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d), modelName = "VII"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  sigmasq <- parameters$variance$sigmasq
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  diag(rep(sqrt(sigmasq[k]), d)), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VII")
}

meV <- function(data, z, prior = NULL, control = emControl(), 
                Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal length of data")
  K <- dimz[2]
  if(!is.null(Vinv)) 
    { G <- K - 1
      if (Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) 
    { WARNING <- "z is missing"
      if(warn) warning(WARNING)
      variance <- list(modelName = "V", d=1, G=G, sigmasq = rep(NA,G))
      parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                         variance=variance, Vinv=Vinv)
      return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                            z=z, parameters=parameters, 
                            control=control, loglik=NA), 
                       WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("me1v",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(G),
                     as.double(if(is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(G),
                     double(G),
                     double(K),
                     PACKAGE = "mclust")[6:12]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "V"), prior[names(prior) !=
                                                                                               "functionName"]))
    temp <- .Fortran("me1vp",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(G),
                     double(G),
                     double(K),
                     PACKAGE = "mclust")[c(10:16, 9)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- temp[[5]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[6]]
  pro <- temp[[7]]
  ## logpost <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6) || 
     any(sigmasq <= max(control$eps, 0))) 
    { WARNING <- "sigma-squared falls below threshold"
      if(warn) warning(WARNING)
      mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
      ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else if(its >= control$itmax[1]) {
    WARNING <- "iteration limit reached"
    if(warn) warning(WARNING)
    its <-  - its
    ret <- 1
  }
  else ret <- 0
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(names(data),NULL)
  variance = list(modelName = "V", d = 1, G = G, 
                  sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik),
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one-dimensional")
  data <- as.vector(data)
  n <- length(data)
  ##
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  # number of groups 
  G <- dimz[2]
  ##
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "V", d=1, G=G, sigmasq=rep(NA,G))
    parameters <- list(pro=rep(NA,G), mean=rep(NA,G), 
                       variance=variance)
    return(structure(list(modelName="V", prior=prior, n=n, d=1, G=G,
                          z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("ms1v",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(G),
                     double(G),
                     double(G),
                     double(G),
                     PACKAGE = "mclust")[5:7]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "V"), prior[names(prior) !=
                                                                                               "functionName"]))
    storage.mode(z) <- "double"
    temp <- .Fortran("ms1vp",
                     as.double(data),
                     z,
                     as.integer(n),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(G),
                     double(G),
                     double(G),
                     PACKAGE = "mclust")[9:11]
  }
  mu <- temp[[1]]
  names(mu) <- as.character(1:G)
  sigmasq <- temp[[2]]
  pro <- temp[[3]]
  WARNING <- NULL
  if(any(sigmasq > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- sigmasq[] <- z[] <- loglik <- NA
    print(G)
    print(sigmasq)
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(names(data),NULL)
  variance = list(modelName = "V", d = 1, G = G, 
                  sigmasq = sigmasq, scale = sigmasq)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "V", prior = prior, n = n, d = 1, G = G, 
                 z = z, parameters = parameters),
            WARNING = WARNING, returnCode = ret)
}

simV <- function(parameters, n, seed = NULL, ...)
{
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, 2), modelName = "V"))
  }
  if(!is.null(seed))
    set.seed(seed)
  mu <- parameters$mean
  G <- length(mu)
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- rep(0, n)
  sd <- sqrt(parameters$variance$sigmasq)
  for(k in 1:G) {
    x[clabels == k] <- mu[k] + rnorm(ctabel[k], sd = sd[k])
  }
  structure(cbind(group = clabels, "1" = x), modelName = "V")
}

cdensVVI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
       any(is.null(parameters[c("pro", "mu", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(z, logarithm = logarithm, modelName = "VVI", 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvvi",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, G)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)  
  structure(z, logarithm = logarithm, modelName = "VVI",
            WARNING = WARNING, retrinCode = ret)
}

emVVI <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepVVI(data, parameters = parameters, warn = warn)$z  
  meVVI(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepVVI <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if (is.null(Vinv) || Vinv <= 0) 
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mu", "variance")]))) ||
       any(is.null(parameters[c("pro", "mu", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VVI", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$scale) ||
        is.null(parameters$variance$shape)) 
    stop("variance parameters are missing")
  temp <- .Fortran("esvvi",
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$scale),
                   as.double(parameters$variance$shape),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:11]
  loglik <- temp[[1]]
  z <- matrix(temp[[2]], n, K)
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL) 
  structure(list(modelName = "VVI", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meVVI <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) > 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv)) {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VVI", d = p, G = G, 
                     scale = rep(NA,G), shape = matrix(as.double(NA),p,G)) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="VVI", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvi",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[7:14]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VVI"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvip",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(K),
                     PACKAGE = "mclust")[11:18]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  scale <- temp[[6]]
  shape <- matrix(temp[[7]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[8]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    sigma <- array(NA, c(p, p, G))
    mu[] <- pro[] <- z[] <- loglik <- shape[] <- NA
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = err)
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VVI", d = p, G = G, 
                   sigma = sigma, scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv)
  structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepVVI <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VII", d=p, G=G, sigmasq=rep(NA,G))
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="VII", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters), 
                     WARNING = WARNING, returnCode = 9))
    
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("msvvi",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(G),
                     PACKAGE = "mclust")[6:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VVI"), prior[names(
                                                          prior) != "functionName"]))
    temp <- .Fortran("msvvip",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(p * G),
                     double(G),
                     double(p * G),
                     double(G),
                     PACKAGE = "mclust")[10:13]
  }
  mu <- matrix(temp[[1]], p, G)
  scale <- temp[[2]]
  shape <- matrix(temp[[3]], p, G)
  dimnames(mu) <- dimnames(shape) <- list(NULL, as.character(1:G))
  pro <- temp[[4]]
  WARNING <- NULL
  if(any(c(scale, shape) > signif(.Machine$double.xmax, 6)) || any(!
                                                                     c(scale, shape))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- shape <- scale[] <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else {
    sigma <- array(apply(sweep(shape, MARGIN = 2, STATS = scale,
                               FUN = "*"), 2, diag), c(p, p, G))
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  variance <- list(modelName = "VVI", d = p, G = G, 
                   sigma = sigma, sigmasq = scale, 
                   scale = scale, shape = shape)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VVI", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters),
            WARNING = WARNING, returnCode = ret)
}

simVVI <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "VVI"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  rtshape <- sqrt(parameters$variance$shape)
  if(!all(dim(rtshape) == dim(mu)))
    stop("shape incompatible with mean")
  rtscale <- sqrt(parameters$variance$scale)
  if(length(rtscale) != G)
    stop("scale incompatible with mean")
  for(k in 1:G) {
    m <- ctabel[k]
    x[clabels == k,  ] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% 
                                  diag(rtscale[k] * rtshape[, k]), MARGIN = 2, STATS = mu[, k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VVI")
}

cdensVVV <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,G)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$cholsigma))
    stop("variance parameters are missing")
  temp <- .Fortran("esvvv",
                   as.logical(1),
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$cholsigma),
                   as.double(-1),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(-1),
                   double(p),
                   double(1),
                   double(n * G),
                   PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, G)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      if(warn) warning(WARNING)
    }
    z[] <- NA
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- NA
    ret <- -1
  }
  else {
    if (!logarithm) z <- exp(z)
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(z, logarithm = logarithm, modelName = "VVV",
            WARNING = WARNING, returnCode = ret)
}

emVVV <- function(data, parameters, prior = NULL, control = emControl(), 
                  warn = NULL, ...)
{
  z <- estepVVV(data, parameters = parameters, warn = warn)$z  
  meVVV(data, z = z, prior = prior, control = control, 
        Vinv = parameters$Vinv, warn = warn)
}

estepVVV <- function(data, parameters, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  if(is.null(dimdat) || length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  pro <- parameters$pro
  pro <- pro/sum(pro)
  l <- length(pro)
  mu <- as.matrix(parameters$mean)
  G <- ncol(mu)
  noise <- l == G + 1
  if(!noise) {
    if(l != G)
      stop("pro improperly specified")
    K <- G
    Vinv <- NULL
  }
  else {
    K <- G + 1
    Vinv <- parameters$Vinv
    if(is.null(Vinv) || Vinv <= 0)
      Vinv <- hypvol(data, reciprocal = TRUE)
  }
  if(any(is.na(unlist(parameters[c("pro", "mean", "variance")]))) ||
       any(is.null(parameters[c("pro", "mean", "variance")]))) {
    WARNING <- "parameters are missing"
    if(warn) warning(WARNING)
    z <- matrix(as.double(NA),n,K)
    dimnames(z) <- list(dimnames(data)[[1]], NULL)
    return(structure(list(modelName = "VVV", n=n, d=p, G=G, z=z,
                          parameters=parameters, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if (is.null(parameters$variance$cholsigma))
    stop("variance parameters are missing")
  temp <- .Fortran("esvvv",
                   as.logical(1),
                   as.double(data),
                   as.double(mu),
                   as.double(parameters$variance$cholsigma),
                   as.double(pro),
                   as.integer(n),
                   as.integer(p),
                   as.integer(G),
                   as.double(if (is.null(Vinv)) -1 else Vinv),
                   double(p),
                   double(1),
                   double(n * K),
                   PACKAGE = "mclust")[10:12]
  lapackCholInfo <- temp[[1]][1]
  loglik <- temp[[2]]
  z <- matrix(temp[[3]], n, K)
  WARNING <- NULL
  if(lapackCholInfo) {
    if(lapackCholInfo > 0) {
      WARNING <- "sigma is not positive definite"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "input error for LAPACK DPOTRF"
      if(warn) warning(WARNING)
    }
    z[] <- loglik <- NA
    ret <- -9
  }
  else if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "cannot compute E-step"
    if(warn) warning(WARNING)
    z[] <- loglik <- NA
    ret <- -1
  }
  else ret <- 0
  dimnames(z) <- list(dimnames(data)[[1]],NULL)
  structure(list(modelName = "VVV", n = n, d = p, G = G, 
                 z = z, parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret)
}

meVVV <- function(data, z, prior = NULL, control = emControl(), 
                  Vinv = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should in the form of a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("data and z should have the same row dimension")
  K <- dimz[2]
  if (!is.null(Vinv))  {
    G <- K - 1
    if(Vinv <= 0) Vinv <- hypvol(data, reciprocal = TRUE)
  }
  else G <- K
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VVV", d = p, G = G, 
                     sigma = array(NA, c(p,p,G)), cholsigma = array(NA, c(p,p,G))) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance, Vinv=Vinv)
    return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          control=control, loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  storage.mode(z) <- "double"
  if(is.null(prior)) {
    temp <- .Fortran("mevvv",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(p * p * G),
                     double(K),
                     double(p),
                     double(p*p),
                     PACKAGE = "mclust")[7:13]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VVV"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("mevvvp",
                     as.logical(control$equalPro),
                     as.double(data),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(if (is.null(Vinv)) -1 else Vinv),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     z,
                     as.integer(control$itmax[1]),
                     as.double(control$tol[1]),
                     as.double(control$eps),
                     double(p * G),
                     double(p * p * G),
                     double(K),
                     double(p),
                     double(p*p),
                     PACKAGE = "mclust")[c(11:17, 10)]
  }
  z <- temp[[1]]
  its <- temp[[2]]
  err <- temp[[3]]
  loglik <- temp[[4]]
  mu <- matrix(temp[[5]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- array(temp[[6]], c(p, p, G))
  pro <- temp[[7]]
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "singular covariance"
    if(warn) warning(WARNING)
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- -1
  }
  else if(loglik <  - signif(.Machine$double.xmax, 6)) {
    if(control$equalPro) {
      WARNING <- "z column sum fell below threshold"
      if(warn) warning(WARNING)
    }
    else {
      WARNING <- "mixing proportion fell below threshold"
      if(warn) warning(WARNING)
    }
    mu[] <- pro[] <- z[] <- loglik <- NA
    sigma <- array(NA, c(p, p, G))
    ret <- if(control$equalPro) -2 else -3
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                   c(p,p,G))
    if(its >= control$itmax[1]) {
      WARNING <- "iteration limit reached"
      if(warn) warning(WARNING)
      its <-  - its
      ret <- 1
    }
    else ret <- 0
  }
  info <- c(iterations = its, error = abs(err))
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- dimnames(cholsigma) <-
    list(dimnames(data)[[2]], dimnames(data)[[2]], NULL)
  variance <- list(modelName = "VVV", d = p, G = G,
                   sigma = sigma, cholsigma = cholsigma)
  parameters <- list(pro=pro, mean=mu, variance=variance, Vinv=Vinv) 
  structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters, control = control,
                 loglik = loglik), 
            info = info, WARNING = WARNING, returnCode = ret)
}

mstepVVV <- function(data, z, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD || length(dimdat) != 2)
    stop("data should be a matrix or a vector")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  z <- as.matrix(z)
  dimz <- dim(z)
  if(dimz[1] != n)
    stop("row dimension of z should equal data length")
  G <- dimz[2]
  if(all(is.na(z))) {
    WARNING <- "z is missing"
    if(warn) warning(WARNING)
    variance <- list(modelName = "VVV", d = p, G = G, 
                     sigma <- array(NA, c(p,p, G)), 
                     cholsigma = array(NA, c(p,p,G))) 
    parameters <- list(pro=rep(NA,G), mean=matrix(as.double(NA),p,G), 
                       variance=variance)
    return(structure(list(modelName="VVV", prior=prior, n=n, d=p, 
                          G=G, z=z, parameters=parameters, 
                          loglik=NA), 
                     WARNING = WARNING, returnCode = 9))
  }
  if(any(is.na(z)) || any(z < 0) || any(z > 1))
    stop("improper specification of z")
  if(is.null(prior)) {
    temp <- .Fortran("msvvv",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     double(p),
                     double(p * G),
                     double(p * p * G),
                     double(G),
                     double(p * p),
                     PACKAGE = "mclust")[7:9]
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = G, modelName = "VVV"), 
                                                 prior[names(prior) != "functionName"]))
    temp <- .Fortran("msvvvp",
                     as.double(data),
                     as.double(z),
                     as.integer(n),
                     as.integer(p),
                     as.integer(G),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(if(any(priorParams$scale != 0)) chol(priorParams$
                                                                      scale) else priorParams$scale),
                     as.double(priorParams$dof),
                     double(p),
                     double(p * G),
                     double(p * p * G),
                     double(G),
                     double(p * p),
                     PACKAGE = "mclust")[11:13]
  }
  mu <- matrix(temp[[1]], p, G)
  dimnames(mu) <- list(NULL, as.character(1:G))
  cholsigma <- array(temp[[2]], c(p, p, G))
  pro <- temp[[3]]
  WARNING <- NULL
  if(any(c(mu, cholsigma) > signif(.Machine$double.xmax, 6))) {
    WARNING <- "cannot compute M-step"
    if(warn) warning(WARNING)
    mu[] <- sigma[] <- cholsigma[] <- NA
    ret <- -1
  }
  else {
    sigma <- array(apply(cholsigma, 3, unchol, upper = TRUE), 
                   c(p,p,G))
    ret <- 0
  }
  dimnames(z) <- list(dimnames(data)[[1]], NULL)
  dimnames(mu) <- list(dimnames(data)[[2]], NULL)
  dimnames(sigma) <- list(dimnames(data)[[2]], dimnames(data)[[2]],
                          NULL)
  dimnames(cholsigma) <- list(dimnames(data)[[2]], 
                              dimnames(data)[[2]], NULL)
  variance <- list(modelName = "VVV", d = p, G = G, 
                   sigma = sigma, cholsigma= cholsigma)
  parameters <- list(pro=pro, mean=mu, variance=variance)
  structure(list(modelName = "VVV", prior = prior, n = n, d = p, G = G, 
                 z = z, parameters = parameters), 
            WARNING = WARNING, returnCode = ret)
  
}

simVVV <- function(parameters, n, seed = NULL, ...)
{
  if(!is.null(seed)) set.seed(seed)
  mu <- as.matrix(parameters$mean)
  d <- nrow(mu)
  G <- ncol(mu)
  if(any(is.na(parameters[c("mean", "variance")])) || 
     any(is.null(parameters[c("mean", "variance")]))) 
    { warning("parameters are missing")
      return(structure(matrix(as.double(NA), n, d + 1), modelName = "VVV"))
  }
  pro <- parameters$pro
  if(is.null(pro))
    pro <- rep(1/G, G)
  clabels <- sample(1:G, size = n, replace = TRUE, prob = pro)
  ctabel <- tabulate(clabels, nbins=G)
  x <- matrix(0, n, d)
  if(is.null(cholsigma <- parameters$variance$cholsigma)) {
    if(is.null(sigma <- parameters$variance$sigma)) {
      stop("variance parameters must inlcude either sigma or cholsigma"
      )
    }
    cholsigma <- apply(sigma, 3, chol)
    for(k in 1:ncol(cholsigma))
      sigma[,  , k] <- cholsigma[, k]
    cholsigma <- sigma
  }
  if(dim(cholsigma)[3] != G)
    stop("variance incompatible with mean")
  for(k in 1:G) 
     { m <- ctabel[k]
       x[clabels == k,] <- sweep(matrix(rnorm(m * d), nrow = m, ncol = d) %*% cholsigma[,,k], 
                                 MARGIN = 2, STATS = mu[,k], FUN = "+")
  }
  dimnames(x) <- list(NULL, paste0("x", 1:d))
  structure(cbind(group = clabels, x), modelName = "VVV")
}


# single component univariate case
mvnX <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(!oneD)
    stop("data must be one dimensional")
  data <- as.vector(data)
  n <- length(data)
  if(is.null(prior)) {
    temp <- .Fortran("mvn1d",
                     as.double(data),
                     as.integer(n),
                     double(1),
                     double(1),
                     double(1),
                     PACKAGE = "mclust")[3:5]
    logpost <- NULL
  }
  else {
    priorParams <- do.call(prior$functionName, c(list(data = 
                                                        data, G = 1, modelName = "X"), prior[names(prior) !=
                                                                                               "functionName"]))
    temp <- .Fortran("mvn1p",
                     as.double(data),
                     as.integer(n),
                     as.double(priorParams$shrinkage),
                     as.double(priorParams$mean),
                     as.double(priorParams$scale),
                     as.double(priorParams$dof),
                     double(1),
                     double(1),
                     double(1),
                     PACKAGE = "mclust")[c(7:9, 6)]
    logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  sigmasq <- temp[[2]]
  loglik <- temp[[3]]
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) {
    WARNING <- "sigma-squared vanishes"
    if(warn) warning(WARNING)
    loglik <- NA
    ret <- -1
  }
  variance = list(modelName= "X", d = 1, G = 1, sigmasq = sigmasq)
  parameters <- list(pro = 1, mean = mu, variance = variance)
  structure(list(modelName = "X", prior = prior, n = n, d = 1, G = 1, 
                 parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret) 
}

cdensX <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- as.name("cdensE")
  z <- eval(mc, parent.frame())
  attr(z, "modelName") <- "X"
  return(z)
}

emX <- function(data, prior = NULL, warn = NULL, ...)
{
  mvnX(data, prior = prior, warn = warn, ...)
}

meX <- emX

# single component multivariate case with diagonal common variance
mvnXII <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) 
    { temp <- .Fortran("mvnxii",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       double(p),
                       double(1),
                       double(1),
                       PACKAGE = "mclust")[4:6]
      logpost <- NULL
  }
  else 
    { priorParams <- do.call(prior$functionName, c(list(data = data, 
                                                        G = 1, 
                                                        modelName = "XII"),
                                                   prior[names(prior) != "functionName"]))
      temp <- .Fortran("mnxiip",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       as.double(priorParams$shrinkage),
                       as.double(priorParams$mean),
                       as.double(priorParams$scale),
                       as.double(priorParams$dof),
                       double(p),
                       double(1),
                       double(1),
                       PACKAGE = "mclust")[c(8:10, 7)]
      logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  sigmasq <- temp[[2]]
  loglik <- temp[[3]]
  Sigma <- sigmasq * diag(p)
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) 
    { WARNING <- "singular covariance"
      if(warn) warning(WARNING)
      loglik <- NA
      ret <- -1
  }
  variance <- list(modelName = "XII", d = p, G = 1,  
                   sigmasq   = sigmasq, Sigma = Sigma, 
                   sigma = array(Sigma, c(p, p, 1)), scale = sigmasq)
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), variance = variance) 
  structure(list(modelName = "XII", prior = prior, n = n, d = p, G = 1, 
                 parameters = parameters, loglik = loglik), 
            WARNING = WARNING, returnCode = ret) 
}

cdensXII <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- as.name("cdensEII")
  z <- eval(mc, parent.frame())
  attr(z, "modelName") <- "XII"
  return(z)
}

emXII <- function(data, prior = NULL, warn = NULL, ...)
{
  mvnXII(data, prior = prior, warn = warn, ...)
}

meXII <- emXII

# single component multivariate case with diagonal different variances
mvnXXI <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) 
    { temp <- .Fortran("mvnxxi",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       double(p),
                       double(1),
                       double(p),
                       double(1),
                       PACKAGE = "mclust")[4:7]
      logpost <- NULL
  }
  else 
    { priorParams <- do.call(prior$functionName, c(list(data = data, 
                                                        G = 1, 
                                                        modelName = "XXI"),
                                                   prior[names(prior) != "functionName"]))
      temp <- .Fortran("mnxxip",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       as.double(priorParams$shrinkage),
                       as.double(priorParams$mean),
                       as.double(priorParams$scale),
                       as.double(priorParams$dof),
                       double(p),
                       double(1),
                       double(p),
                       double(1),
                       PACKAGE = "mclust")[c(8:11, 7)]
      logpost <- temp[[5]]
  }
  mu <- temp[[1]]
  scale <- temp[[2]]
  shape <- temp[[3]]
  loglik <- temp[[4]]
  Sigma <- diag(scale * shape)
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) 
    { WARNING <- "singular covariance"
      if(warn) warning(WARNING)
      loglik <- NA
      ret <- -1
  }
  variance <- list(modelName = "XXI", d = p, G = 1, 
                   Sigma = Sigma, sigma = array(Sigma, c(p, p, 1)),
                   scale = scale, shape = shape)
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), 
                     variance = variance)
  structure(list(modelName = "XXI", prior = prior, n = n, d = p, G = 1, 
                 parameters = parameters, loglik = loglik),
            WARNING = WARNING, returnCode = ret) 
}

cdensXXI <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- as.name("cdensEEI")
  z <- eval(mc, parent.frame())
  attr(z, "modelName") <- "XXI"
  return(z)
}

emXXI <- function(data, prior = NULL, warn = NULL, ...)
{
  mvnXXI(data, prior = prior, warn = warn, ...)
}

meXXI <- emXXI

# single component multivariate case with full covariance matrix
mvnXXX <- function(data, prior = NULL, warn = NULL, ...)
{
  if(is.null(warn)) warn <- mclust.options("warn")
  dimdat <- dim(data)
  oneD <- is.null(dimdat) || length(dimdat[dimdat > 1]) == 1
  if(oneD)
    stop("for multidimensional data only")
  if(length(dimdat) != 2)
    stop("data must be a matrix")
  data <- as.matrix(data)
  n <- nrow(data)
  p <- ncol(data)
  if(is.null(prior)) 
    { temp <- .Fortran("mvnxxx",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       double(p),
                       double(p * p),
                       double(1),
                       PACKAGE = "mclust")[c(4:6)]
      logpost <- NULL
  }
  else 
    { priorParams <- do.call(prior$functionName, c(list(data = data, 
                                                        G = 1, 
                                                        modelName = "XXX"),
                                                   prior[names(prior) != "functionName"]))
      temp <- .Fortran("mnxxxp",
                       as.double(data),
                       as.integer(n),
                       as.integer(p),
                       double(p),
                       as.double(priorParams$shrinkage),
                       as.double(priorParams$mean),
                       as.double(if(any(priorParams$scale != 0)) 
                                   chol(priorParams$scale) else priorParams$scale),
                       as.double(priorParams$dof),
                       double(p),
                       double(p * p),
                       double(1),
                       PACKAGE = "mclust")[c(9:11, 8)]
      logpost <- temp[[4]]
  }
  mu <- temp[[1]]
  cholSigma <- matrix(temp[[2]], p, p)
  Sigma <- unchol(cholSigma, upper = TRUE)
  loglik <- temp[[3]]
  ## Sigma = t(cholSigma) %*% cholSigma
  ret <- 0
  WARNING <- NULL
  if(loglik > signif(.Machine$double.xmax, 6)) 
    { WARNING <- "singular covariance"
      if(warn) warning(WARNING)
      loglik <- NA
      ret <- -1
  }
  variance <- list(modelName = "XXX", d = p, G = 1,
                   Sigma = Sigma, cholSigma = cholSigma, 
                   cholsigma = cholSigma,
                   sigma = array(Sigma, c(p, p, 1))) 
  parameters <- list(pro = 1, mean = matrix(mu, ncol = 1), 
                     variance = variance)
  structure(list(modelName = "XXX", prior = prior, n = n, d = p, G = 1,
                 parameters = parameters, loglik = loglik), 
            WARNING = WARNING, returnCode = ret)
}

cdensXXX <- function(data, logarithm = FALSE, parameters, warn = NULL, ...)
{
  call <- match.call()
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- as.name("cdensEEE")
  z <- eval(mc, parent.frame())
  attr(z, "modelName") <- "XXX"
  return(z)
}

emXXX <- function(data, prior = NULL, warn = NULL, ...)
{
  mvnXXX(data, prior = prior, warn = warn, ...)
}

meXXX <- emXXX