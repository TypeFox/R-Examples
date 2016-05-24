mlogit <- function(formula, data, subset, weights, na.action, start= NULL,
                   alt.subset = NULL, reflevel= NULL,
                   nests = NULL, un.nest.el = FALSE, unscaled = FALSE,
                   heterosc = FALSE, rpar = NULL, probit = FALSE,
                   R = 40, correlation = FALSE, halton = NULL, random.nb = NULL,
                   panel = FALSE, estimate = TRUE, seed = 10, ...){

  start.time <- proc.time()
  callT <- match.call(expand.dots = TRUE)
  callF <- match.call(expand.dots = FALSE)
  formula <- callF$formula <- mFormula(formula)
  nframe <- length(sys.calls())
  heterosc.logit <- heterosc
  nested.logit <- !is.null(nests)
  if (!is.null(nests) && length(nests) == 1 && nests == "pcl"){
    nested.logit <- FALSE
    pair.comb.logit <- TRUE
  }
  else pair.comb.logit <- FALSE
  if (nested.logit) attr(nests, "unique") <- ifelse(un.nest.el, TRUE, FALSE)
  mixed.logit <- !is.null(rpar)
  multinom.logit <- !heterosc & is.null(nests) & is.null(rpar) & !probit
  if (heterosc.logit + nested.logit + mixed.logit + probit > 1)
    stop("only one of heterosc, rpar, nests and probit can be used")
  if (multinom.logit){
    if (is.null(callT$method)) callT$method <- 'nr'
    if (is.null(callT$print.level)) callT$print.level <- 0
  }
  
  # 1 ############################################################
  #  check whether arguments for mlogit.data are present: if so run
  #  mlogit.data
  ################################################################

  mldata <- callT
  response.name <- paste(deparse(attr(formula, "lhs")[[1]]))
  m <- match(c("data", "choice", "shape", "varying", "sep",
               "alt.var", "chid.var", "alt.levels",
               "opposite", "drop.index", "id", "ranked"),
             names(mldata), 0L)
  use.mlogit.data <- sum(m[-1]) > 0
  if (use.mlogit.data){
    mldata <- mldata[c(1L, m)]
    mldata[[1L]] <- as.name("mlogit.data")
    mldata$choice <- response.name
#    data <- eval(mldata, sys.frame(which = nframe))
    data <- eval(mldata, parent.frame())
  }

  # 2 ###########################################################
  # model.frame (with the data provided or the one computed by
  # mlogit.data
  ###############################################################

  mf <- callT
  m <- match(c("formula", "data", "subset", "na.action", "weights"),
             names(mf), 0L)
  mf <- mf[c(1L, m)]
  # Est-ce important ??
#  mf$drop.unused.levels <- TRUE
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf$data <- data # fix the bug when the data is called mldata
#  mf <- eval(mf, sys.frame(which = nframe))
  mf <- eval(mf, parent.frame())

  # change the reference level of the response if required
  if (!is.null(reflevel)) attr(mf, "index")[["alt"]] <- relevel(attr(mf, "index")[["alt"]], reflevel)
  index <- attr(mf, "index")
  alt <- index[["alt"]]
  chid <- index[["chid"]]
  if (panel){
    if (!mixed.logit) stop("panel is only relevant for mixed logit models")
    id <- index[["id"]]
    if (is.null(id)) stop("no individual index")
    id <- split(index[["id"]], alt)[[1]]
  }
  else id <- NULL
  # compute the relevent subset if required
  if (!is.null(alt.subset)){
    # we keep only choices that belong to the subset
    choice <- alt[model.response(mf)]
    choice <- choice %in% alt.subset
    unid <- unique(chid)
    names(choice) <- as.character(unid)
    id.kept <- choice[as.character(chid)]
    # we keep only the relevant alternatives
    alt.kept <- alt %in% alt.subset
    # the relevant subset for the data.frame and the indexes
    mf <- mf[id.kept & alt.kept, , drop = FALSE]
    alt <- alt[id.kept & alt.kept , drop = TRUE]
    chid <- chid[id.kept & alt.kept , drop = TRUE]
  }
  # balanced the data.frame i.e. insert rows with NA when an
  # alternative is not relevant
  alt.un <- unique(alt)
  chid.un <- unique(chid)
  alt.lev <- levels(alt)
  J <- length(alt.lev)
  n <- length(chid.un)
  T <- length(alt.un)
  balanced <- TRUE
  if (nrow(mf) != (n * T)){
    # la ligne ci-dessous ne sert à rien ?
    rownames(mf) <- paste(chid, alt, sep = ".")
    all.rn <- as.character(t(outer(chid.un, alt.un, paste, sep = ".")))
    mf <- mf[all.rn, ]
    rownames(mf) <- all.rn
    chid <- rep(chid.un, each = T)
    alt <- rep(alt.un, n)
    index <- data.frame(chid = chid, alt = alt, row.names = rownames(mf))
    balanced <- FALSE
  }
  #suppress individuals for which no choice is made
  if (FALSE){
    delete.id <- tapply(model.response(mf), chid, sum, na.rm = TRUE)
    delete.id <- names(delete.id[delete.id == 0])
    mf <- mf[!(chid %in% delete.id), ]
    index <- index[rownames(mf),]
    index[[1]] <- index[[1]][, drop=TRUE]
    index[[2]] <- alt <- index[[2]][, drop=TRUE]
    attr(mf, "index") <- index
  }
  # if estimate is FALSE, return the data.frame
  if (!estimate) return(mf)

  # 3 ###########################################################
  # extract the elements of the model
  ###############################################################

  y <- model.response(mf)
  choice <- na.omit(alt[y])
  X <- model.matrix(formula, mf)
  K <- ncol(X)
  df.residual <- n - K
  colnamesX <- colnames(X)
  if (any(names(mf)=="(weights)")){
    weights <- mf[["(weights)"]] <- mf[["(weights)"]]/mean(mf[["(weights)"]])
    weights <- split(weights, alt)[[1]]
  }
  else weights <- NULL
  freq <- table(alt[y])
  # Xl and yl are lists of length J which contains n matrix / vector
  # of covariates and response (a boolean) ; yv is a vector that
  # contains the chosen alternative
  Xl <- vector(length = J, mode = "list")
  names(Xl) <- levels(alt)
  for (i in levels(alt)){
    Xl[[i]] <- X[alt == i, , drop = FALSE]
  }
  yl <- split(y, alt)
  yl <- lapply(yl, function(x){x[is.na(x)] <- FALSE ; x})

  
  if (probit){
    # for probit the response is a vector that contains the chosen
    # alternative as a numeric ; NA values of y are replaced by FALSE
    # so that the chosen alternative is correctly returned
    y2 <- y
    y2[is.na(y2)] <- FALSE
    yv <- as.numeric(alt[y2])
    # DX is a list of covariates first differenced with respect with the
    # chosen alternative
    DX <- vector("list", length = J - 1)
    DX <- lapply(DX, function(x) return(matrix(NA, nrow = n, ncol = K)))
    
    for (i in 1:n){
      any <- yv[i]
      j <- 1
      newj <- 1
      for (j in 1:J){
        if (j != any){
          DX[[newj]][i,] <- Xl[[j]][i, ] - Xl[[any]][i,]
          newj <- newj + 1
        }
      }
    }
  }

  # 4 ###########################################################
  # Compute the starting values
  ###############################################################

  # first give names to the supplementary coefficients and values if
  # start is null
  if (nested.logit){
    L <- length(nests)
    if (un.nest.el){
      if (is.null(start) || length(start) == K) sup.coef <- c(iv = 1)
      names.sup.coef <- 'iv'
    }
    else{
      if (is.null(start) || length(start) == K) sup.coef <- rep(1, L)
      names.sup.coef <- paste("iv", names(nests), sep = ".")
    }
  }
  if (pair.comb.logit){
    if (un.nest.el){
      if (is.null(start)) sup.coef <- c(iv = 1)
      names.sup.coef <- 'iv'
    }
    else{
      names.sup.coef <- NULL
      for (i in 1:(J-1)){
        names.sup.coef <- c(names.sup.coef, paste('iv', alt.lev[i], alt.lev[(i+1):J], sep = "."))
      }
      sup.coef <- rep(1, length(names.sup.coef))
      # in case we have to suppress one nest
      ## sup.coef <- sup.coef[-1]
      ## names.sup.coef <- names.sup.coef[-1]
    }
  }
  if (heterosc.logit){
    if (is.null(start) || length(start) == K) sup.coef <- rep(1, J-1)
    names.sup.coef <- paste("sp", alt.lev[-1], sep = ".")
  }
  if (mixed.logit){
      unknowndist <- rpar[! (rpar %in% c("cn", "ln", "tn", "n", "u", "t"))]
      if (length(unknowndist) > 0){
          udstr <- paste("unknown distribution", paste(unique(unknowndist), collapse = ", "))
          stop(udstr)
      }
      Vara <- sort(match(names(rpar), colnamesX))
      Varc <- (1:K)[- Vara]
      nmean <- length(c(Varc, Vara))
      nvar <- length(Vara)
      if (!correlation){
          if (is.null(start) || length(start) == K) sup.coef <- rep(.1, nvar)
          names.sup.coef <- paste("sd", colnamesX[Vara], sep = ".")
      }
      else{
          if (is.null(start) || length(start) == K) sup.coef <- rep(.1, 0.5 * nvar * (nvar + 1))
          names.sup.coef <- c()
          Ka <- length(rpar)
          for (i in 1:Ka){
              names.sup.coef <- c(names.sup.coef,
                                  paste(names(rpar)[i], names(rpar)[i:Ka], sep = "."))
          }
      }
  }
  if (probit){
    names.sup.coef <- c()
    for (i in 2:J){
      names.sup.coef <- c(names.sup.coef,
                          paste(alt.lev[i], alt.lev[i:J], sep = "."))
    }
    if (is.null(start) || length(start) == K){
      corrMat <- matrix(0.5, J - 1, J - 1)
      diag(corrMat) <- 1
      corrMat <- chol(corrMat)
      sup.coef <- (t(corrMat)[!upper.tri(corrMat)])
      names(sup.coef) <- names.sup.coef
    }
  }
  ## start can be :
  ##   1. NULL, in this case estimate the multinomial logit model,
  ##   2. a vector of length K ; then add starting values for the
  ##   supplementary coefficients,
  ##   3. a full set ; then just name the coefs.
  if (is.null(start)){
    callst <- callT
    start <- rep(0, K)
    names(start) <- colnamesX
    callst$start <- start
    callst$print.level <- 0
    if (!multinom.logit){
      callst$nests <- callst$rpar <- callst$constPar <- callst$iterlim <- NULL
      callst$heterosc <- callst$panel <- callst$probit <- FALSE
      callst$print.level <- 0
#      start <- coef(eval(callst, sys.frame(which=nframe)))
      start <- coef(eval(callst, parent.frame()))
      if (mixed.logit){
        ln <- names(rpar[rpar == "ln"])
        start[ln] <- log(start[ln])
      }
    }
  }
  if (length(start) == K){
    names(start) <- colnamesX
    if (!multinom.logit){
      names(sup.coef) <- names.sup.coef
      if (probit) start <- start * sqrt(3) / pi
      start <- c(start, sup.coef)
    }
  }
  else{
    if (!multinom.logit) names(start) <- c(colnamesX, names.sup.coef)
    else names(start) <- colnamesX
  }

  # if constPar is numeric, insert the relevant value in the start
  # vector and transform the constPar vector to character

  # 5 ###################################################################
  # Estimate the model using mlogit.nlm and passing the correct arguments
  #######################################################################

  opt <- callT
  if (!is.null(opt$constPar)){
    theconstPar <- eval(opt$constPar)
    if (is.numeric(theconstPar)){
      if (is.null(names(theconstPar))) stop('the numeric constPar vector should be named')
      start[names(theconstPar)] <- theconstPar
      opt$constPar <- names(theconstPar)
    }
  }
  if (probit) if (is.null(opt$constPar)) opt$constPar <- names.sup.coef[1]
  opt$start <- start
  m <- match(c("method", "print.level", "iterlim",
               "start", "constPar","tol", "ftol", "steptol"),
             names(opt), 0L)
  opt <- opt[c(1L, m)]
  opt[[1]] <- as.name('mlogit.optim')
  opt$logLik <- as.name('lnl.slogit')

  opposite <- TRUE
  if (is.null(weights)) weights <- 1
  opposite <- ifelse(opposite, -1, +1)
  opt[c('weights', 'opposite')] <- list(as.name('weights'), as.name('opposite'))
  if (!probit) opt[c('X', 'y')] <- list(as.name('Xl'), as.name('yl'))
  else opt[c('X', 'y')] <- list(as.name('DX'), as.name('yv'))
  if (mixed.logit){
    opt$logLik <- as.name('lnl.rlogit')
    opt[c('R', 'seed', 'id', 'rpar', 'correlation', 'halton')] <-
      list(as.name('R'), as.name('seed'), as.name('id'), as.name('rpar'), as.name('correlation'), as.name('halton'))
  }
  if (probit){
    opt$logLik <- as.name('lnl.mprobit')
    opt[c('R', 'seed')] <-
      list(as.name('R'), as.name('seed'))
  }
  if (heterosc.logit){
    opt$logLik <- as.name('lnl.hlogit')
    rn <- gauss.quad(R, kind = "laguerre")
    opt[c('rn')] <- list(as.name('rn'))
  }
  if (nested.logit){
    opt$logLik <- as.name('lnl.nlogit')
    opt$nests <- as.name('nests')
    opt$un.nest.el <- as.name('un.nest.el')
    opt$unscaled <- as.name('unscaled')
  }
  if (pair.comb.logit){
    opt$logLik <- as.name('lnl.nlogit')
    alt1 <- rep(alt.lev, c((J-1):0))
    alt2 <- alt.lev[unlist(lapply(2:J, function(x) x:J))]
    names.nests <- paste(alt1, alt2, sep = ".")
    nests <- mapply(function(x,y) c(x,y), alt1, alt2, SIMPLIFY = FALSE)
    names(nests) <- names.nests
    # in case we have to suppress one nest
    # lnests <- lnests[-1]
    opt$nests <- nests
    opt$unscaled <- as.name('unscaled')
    opt$un.nest.el <- as.name('un.nest.el')
  }
  ## Plante avec parent.frame
  x <- eval(opt, sys.frame(which = nframe))

  # 6 ###########################################################
  # put the result in form
  ###############################################################

  # some general features
  n <- sum(freq)
  x$est.stat$elaps.time <- proc.time() - start.time
  logLik <- structure( - as.numeric(x$optimum),
                      df = length(x$coefficients),
                      null = sum(freq * log(freq / n)),
                      class = "logLik"
                      )
  if (mixed.logit) rpar <- make.rpar(rpar, correlation, x$coefficients, NULL) else rpar <- NULL
  if (!(nested.logit | pair.comb.logit)) nests <- NULL
  
  # if no hessian is returned, use the BHHH approximation
  if (is.null(attr(x$optimum, 'hessian'))) hessian <- - crossprod(attr(x$optimum, 'gradi'))
  else hessian <- - attr(x$optimum, 'hessian')
  fitted <- attr(x$optimum, "fitted")
  probabilities <- attr(x$optimum, "probabilities")
  resid <- Reduce("cbind", yl) - fitted
  attr(x$coefficients, "fixed") <- attr(x$optimum, "fixed")
  gradient <- -  attr(x$optimum, "gradient")
  gradient <- - attr(x$optimum, "gradi")

  
  # compute the probabilities for all the alternatives for
  # heteroscedastic and the probit model
  if (probit | heterosc){
    opt$logLik <- opt$iterlim <- opt$method <- opt$print.level <- opt$tol <- opt$constPar <- NULL
    names(opt)[[2]] <- 'param'
    opt[[2]] <- x$coefficients
    opt$gradient <- FALSE
    opt[[1]] <- as.name('lnl.hlogit')
    if (probit) opt[[1]] <- as.name('lnl.mprobit') else opt[[1]] <- as.name('lnl.hlogit')
    probabilities <- c()
    for (k in 1:J){
      if (probit){
        they <- rep(k, n)
        opt$y <- they
        DX <- vector("list", length = J - 1)
        DX <- lapply(DX, function(x) return(matrix(NA, nrow = n, ncol = K)))
        for (i in 1:n){
          any <- k
          j <- 1
          newj <- 1
          for (j in 1:J){
            if (j != any){
              DX[[newj]][i, ] <- Xl[[j]][i, ] - Xl[[any]][i, ]
              newj <- newj + 1
            }
          }
        }
        opt$X <- DX
      }
      if (heterosc){
        they <- vector(mode='list', length= J)
        they <- lapply(they, function(x) rep(FALSE, n))
        they[[k]] <- rep(TRUE, n)
        names(they) <- names(yl)
        opt$y <- they
      }
      # Plante avec parent.frame
      probabilities <- cbind(probabilities,
                             attr(eval(opt, sys.frame(which = nframe)), 'fitted'))
    }
    colnames(probabilities) <- alt.lev
    attr(x$optimum, "probabilities") <- probabilities
  }

  # Compute the covariance matrix of the errors
  if (multinom.logit | mixed.logit){
    alt.names <- colnames(probabilities)
    J <- length(alt.names)
    Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
    diag(Omega) <- pi^2 / 6
  }
  if (probit){
    S <- matrix(0, J - 1, J - 1)
    S[!upper.tri(S)] <- x$coefficients[-c(1:K)]
    Mi <- list()
    M <- rbind(0, diag(J - 2))
    Mi[[1]] <- diag(J - 1)
    for (i in 2:(J - 1)){
      Mi[[i]] <- cbind(M[, 0:(i - 2), drop = FALSE], -1, M[, ((i - 1):(J - 2))])
    }
    Mi[[J]] <- cbind(M, -1)
    names(Mi) <- alt.lev
    Omega <- lapply(Mi, function(x) tcrossprod(x %*% S))
    for (i in 1:J){
      colnames(Omega[[i]]) <- rownames(Omega[[i]]) <- alt.lev[-i]
    }
  }
  if (nested.logit | pair.comb.logit){
    alt.names <- colnames(probabilities)
    J <- length(alt.names)
    Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
    for (i in 1:length(nests)){
      anEl <- x$coefficients[paste('iv', names(nests)[i], sep = ".")]
      alts <- nests[[i]]
      M <- length(alts)
      for (m in 1:M){
        for (n in 1:M){
          Omega[alts[m], alts[n]] <- (1 - anEl^2)
        }
      }
    }
    diag(Omega) <- 1
    Omega <- Omega * tcrossprod(rep(pi/sqrt(6), J))
  }
  if (heterosc.logit){
    alt.names <- colnames(probabilities)
    J <- length(alt.names)
    Omega <- matrix(0, J, J, dimnames = list(alt.names, alt.names))
    diag(Omega) <- pi^2 / 6
    for (i in 2:J){
      analt <- alt.names[i]
      Omega[analt, analt] <- pi^2 / 6 * x$coefficients[paste('sp', analt, sep = '.')]^2
    }
  }

  result <- structure(
                      list(
                           coefficients  = x$coefficients,
                           logLik        = logLik,
                           gradient      = gradient,
                           hessian       = hessian,
                           est.stat      = x$est.stat,
                           fitted.values = fitted,
                           probabilities = probabilities,
                           residuals     = resid,
                           omega         = Omega,
                           rpar          = rpar,
                           nests         = nests,
                           model         = mf,
                           freq          = freq,
                           formula       = formula,
                           call          = callT),
                      class = 'mlogit'
                      ) 
  #result$Mi <- Mi
  #result$alt.lev <- alt.lev
  result
}
