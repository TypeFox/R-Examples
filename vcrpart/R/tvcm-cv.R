## --------------------------------------------------------- #
##' Author:          Reto Buergin
##' E-Mail:          rbuergin@gmx.ch
##' Date:            2016-02-22
##'
##' Description:
##' Function for model selection and assessment for 'tvcm' objects.
##'
##' Contents:
##' oobloss.tvcm:        computes the out.of-bag loss
##' folds:               parameters for cross-validation folds
##' tvcm_folds:          create cross-validation folds
##' cvloss.tvcm:         cross-validation for 'tvcm' objects
##' print.cvloss.tvcm:   print for 'cv.tvcm' objects
##' plot.cvloss.tvcm:    plot fot 'cv.tvcm' objects
##'
##' Last modifications:
##' 2016-02-22: adapt code of 'tvcm_folds' to allow creating folds
##'             for olmm objects (not officially!!!)
##' 2014-09-09: tvcm_folds: the 'seed' attribute is now the number
##'             of the seed and not the RNG state anymore.
##' 2014-09-07: modifications for direct call from 'tvcm'
##' 2014-09-02: - modifications on 'tvcm_get_node'. The former
##'               implementation was a time-killer an therefore
##'               there is a new argument 'formList' and a
##'               auxiliary function 'tvcm_get_fitted' was added
##' 2014-08-30: - deleted 'nsplit' extension of in 'cvloss'. Reason:
##'               cross-validation for number of split would
##'               require a different pruning procedure
##'             - small justifications for the plot
##' 2014-08-06: - substituted 'cvfolds' by 'folds', which defines
##'               a list of parameters for the new function
##'               'tvcm_folds' that creates the cross-validation
##'               matrix
##'               matrix
##' 2014-07-29: - defined 'cvfolds' as method for 'tvcm' objects
##'             - added new argument 'weights' to 'cvfolds' to
##'               allow for models where the weights represent
##'               counts
##' 2014-07-22: - removed AIC and BIC methods since they do
##'               not apply to the 'tvcm' framework
##' 2014-07-17: - change the column names of the cross-validation
##'               output matrix
##'             - improve the desciption
##' 2014-07-08: - remove the 'sub' argument of AIC.tvcm, BIC.tvcm
##'             - remove additional methods for AIC.tvcm and BIC.tvcm
##'             - remove stabsel.tvcm and and corresponding methods
##'             - cvloss: set 'cp' as the only tuning parameter
##'             - cvloss: add 'direction' as new parameter 
## --------------------------------------------------------- #


oobloss.tvcm <- function(object, newdata = NULL, weights = NULL, 
                         fun = NULL, ...) {
  if (is.null(fun)) {
    fun <- function(y, mu, wt)
      sum(object$info$family$dev.resids(y, mu, wt))
  }
  if (missing(newdata)) stop("require 'newdata'.")
  if (is.null(weights)) weights <- rep(1.0, nrow(newdata))
  yName <- all.vars(object$info$formula$original)[1]
  yMat <- model.matrix(formula(paste("~ -1 + ", yName)), data = newdata)
  if (object$info$family$family == "binomial" && ncol(yMat) > 1L)
    yMat <- yMat[,2L,drop = FALSE]
  mu <- suppressWarnings(predict(object, newdata, type = "response", ...))
  rval <- fun(yMat, mu, weights)
  return(rval)
}


folds_control <- function(type = c("kfold", "subsampling", "bootstrap"),
                          K = ifelse(type == "kfold", 5, 100),
                          prob = 0.5, weights = c("case", "freq"),
                          seed = NULL) {
  if ("bootstrapping" %in% type) type <- "bootstrap"
  type <- match.arg(type)
  stopifnot(is.numeric(K) && length(K) == 1L)
  if (round(K) < 1L) stop("'K' must be a positive number.")
  if (type == "kfold" && K < 2L)
    stop("'K' must be larger than 1 for 'type = 'kfold''")
  if (K != as.integer(round(K)))
    warning(paste("'K' has been set to ", K, ".", sep = ""))
  K <- as.integer(round(K))
  stopifnot(is.numeric(prob))
  if (prob < 0 | prob > 1)
    stop("'prob' must be within the interval [0, 1].")
  K <- round(K)
  stopifnot(is.numeric(prob) && length(prob) == 1L)
  weights <- match.arg(weights)
  return(structure(list(type = type,
                        K = K,
                        prob = prob,
                        weights = weights,
                        seed = seed),
                   class = "folds"))
}


## --------------------------------------------------------- #
##' Creates a cross-validation matrix
##'
##' @param object  an object of class \code{tvcm}
##' @param args    a list of arguments as produced by
##'    \code{\link{folds}}.
##'
##' @return A matrix.
## --------------------------------------------------------- #

tvcm_folds <- function(object, control) {

  stopifnot(inherits(control, "folds"))
  
  ## detach input
  type <- control$type
  K <- control$K
  prob <- control$prob
  weights <- control$weights
  seed <- control$seed

  if (inherits(object, "tvcm")) object <- extract(object, "model")
  
  subject <- object$subject
  if (inherits(object, "olmm")) {
    if (max(table(subject)) > 0L) {
      if (weights == "freq")
        stop("option 'weights = 'freq'' is not available for 'olmm' objects")
      if (type == "bootstrap")
        stop("option 'type = 'bootstrap'' is not available for 'olmm' object",
             "with 2-stage structures.")
    }
  }
  
  freq <- switch(weights,
                 case = rep(1, nobs(object)),
                 freq = weights(object))
  if (weights == "freq" && any(!freq == round(freq)))
    stop("some of the weights are not integers.")
  
  if (is.null(subject)) subject <- factor(rep(1:length(freq), freq))
  N <- nlevels(subject)
    
  getSample <- function(subject, freq, prob, replace, weights) {
    levs <- if (weights == "case") 1:nlevels(subject) else freq = 1:length(subject)
    rSample <- sample(x = levs,
                      size = ceiling(prob * length(levs)),
                      replace = replace)
    rSample <- table(factor(rSample, levels = levs))
    if (weights == "case") rSample <- rSample[subject]
    if (weights == "freq") rSample <- tapply(rSample, subject, sum)
    return(rSample)
  }

  if (!exists(".Random.seed", envir = .GlobalEnv)) runif(1)
  oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
  if (!is.null(seed)) set.seed(seed)
  RNGstate <- .Random.seed
  
  if (type == "subsampling") {
    
    folds <- replicate(K, getSample(subject, freq, prob, FALSE, weights))
    
  } else if (type == "kfold") {
    
    if (type == "kfold" & (K > length(subject)) || (K <= 1)) 
      stop("'K' outside allowable range")

    levs <- if (weights == "case") 1:nlevels(subject) else 1:length(subject)
    selected <- sample(levs, length(levs))
    split <- c(0, quantile(levs, (1:(K - 1)) / K, type = 1L), max(levs))
    folds <- matrix(, length(levs), K)
    for (i in 1:K) {
      subs <- levs > split[i] & levs <= split[i + 1L]
      folds[, i] <- 1.0 * (!levs %in% selected[subs])
    }
    if (weights == "case") folds <- folds[subject, ]
    if (weights == "freq") folds <- apply(folds, 2, tapply, subject, sum)
    
  } else if (type == "bootstrap") {
    
    folds <- replicate(K, getSample(subject, freq, prob = 1, TRUE, weights))
  }

  colnames(folds) <- 1:K
  rownames(folds) <- rownames(model.frame(object))

  assign(".Random.seed", oldSeed, envir=globalenv())
  attr(folds, "type") <- type
  attr(folds, "value") <- ifelse(weights == "freq", "weights", "freq")
  attr(folds, "seed") <- seed
  return(folds)
}


cvloss.tvcm <- function(object, folds = folds_control(), ...) {
  
  mc <- match.call()

  control <- extract(object, "control")
  stopifnot(inherits(folds, "folds"))
  control$folds <- folds
  
  ## set the verbose (no verbose is shown when evaluating the validation sets)
  object$info$control$verbose <- FALSE

  ## desactivate parallelization in single evaluations
  object$info$control$papply <- "lapply"
  object$info$control$papply.args <- list()
  
  ## (hidden) type of evaluation ('loss', 'forest' etc.)
  type <- list(...)$type
  if (is.null(type)) type <- "loss"

  ## (hidden) whether the original sample should be evaluated
  original <- list(...)$original
  if (is.null(original)) original <- FALSE

  ## get folds
  foldsMat <- tvcm_folds(object, folds)

  ## get initial complexity penalty
  cp <- control$cp
  
  ## weights and model frame
  weights <- weights(object$info$model)
  mf <- model.frame(object)
  
  ## cross-validation function

  cvFun <- function(i) {

    ## the return value of 'cvFun'
    cv <- vector(mode = "list", length = switch(type, loss = 2L, forest = 3L))
    
    if (control$verbose) {
      if (control$papply == "lapply" && i > 1L) cat("\n")
      if (control$papply != "lapply") cat("[", i, "]") else cat("* fold", i, "...")
    }
    
    if (i > 0L) {
      ## extract subset
      if (attr(foldsMat, "value") == "weights") {
        ibSubs <- foldsMat[, i] > 0.0
        ibWeights <- foldsMat[ibSubs, i]
        oobSubs <- rep(TRUE, nrow(foldsMat))
        oobWeights <- weights - foldsMat[, i]
      } else {
        ibSubs <- foldsMat[, i]
        ibSubs <- rep(1:length(ibSubs), ibSubs) # rep. bootstrap replications  
        oobSubs <- foldsMat[, i] <= 0
        ibWeights <- weights[ibSubs]
        oobWeights <- weights[oobSubs]
      }
    } else {
      object$info$control <- control
      ibSubs <- NULL
      ibWeights <- NULL
    }
    
    ## re-fit the tree with training sample (in-bag)
    ibTree <- try(tvcm_grow(object, ibSubs, ibWeights), silent = TRUE)
    
    if (i == 0) {
      if (inherits(ibTree, "try-error"))
        stop("partitioning failed.")
      return(ibTree)
    }
    
    if (!inherits(ibTree, "try-error")) {
      
      if (type == "loss") {
        run <- 1L
        
        while (run > 0L) {
          ibTree <- try(prune(tree = ibTree, cp = cp), TRUE)
          if (!inherits(ibTree, "try-error")) {
            ## save the out-of-bag loss and the current tuning parameter
            oobLoss <- oobloss(ibTree, newdata = mf[oobSubs,,drop = FALSE],
                               weights = oobWeights, fun = control$ooblossfun)
            cv[[1L]] <-
              cbind(cv[[1L]], cp)
            cv[[2L]] <-
              cbind(cv[[2L]], oobLoss)

            ## set a new and stronger tuning parameter
            tab <- ibTree$info$prunepath[[length(ibTree$info$prunepath)]]$tab
            if (nrow(tab) > 1L) cp <- min(tab[, "dev"][-1L])
            if (nrow(tab) == 1L) run <- 0L
            
          } else {
            cv <- NULL
            run <- 0L

          }
        }
  
      } else if (type == "forest") {        
        if (control$verbose && control$papply == "lapply") cat("...")
        ibTreePr <- ibTree
        cv[[1]] <- ibTree$info$node
        cv[[2]] <- coef(extract(ibTree, "model"))
        cv[[3]] <- extract(ibTree, "model")$contrasts
        
      }

    } else {
      cv <- ibTree

    }
    if (control$verbose) {
      if (control$papply != "lapply") {
        if (inherits(cv, "try-error")) cat("failed")
      } else {
        if (inherits(cv, "try-error")) cat("failed") else cat(" OK")
      }
    }
    ## return output
    return(cv)
  }

  call <- list(name = as.name(control$papply),
               X = quote(seq(ifelse(original, 0L, 1L), ncol(foldsMat))),
               FUN = quote(cvFun))
  call[names(control$papply.args)] <- control$papply.args
  mode(call) <- "call"
  cv <- eval(call)
  
  if (type %in% c("loss")) {

    ## extract tree on all data
    if (original) {
      tree <- cv[[1L]]
      tree$info$control <- control
      cv <- cv[2:length(cv)]
      if (is.null(tree)) stop("partitioning failed.")
    } else {
      tree <- NULL
    }
    
    ## delete fails
    fails <- sapply(cv, function(x) sum(sapply(x, length)) == 0)
    if (all(fails)) stop("no valid results.")
    cv <- cv[!fails]
    
    if (type == "loss") {

      ## function that evaluates the loss at each alpha of 'grid'
      getVals <- function(x, grid, rowsG = 1L, rowsL = 1L) {
        rval <- matrix(NA, length(row), length(grid))
        for (i in 1:length(grid)) {
          cols <- which(x[[1L]][rowsG,] <= grid[i])
          if (length(cols) > 1)
            cols <- which(x[[1L]][rowsG, ] == max(x[[1L]][rowsG, cols]))
          if (length(cols) > 0)
            rval[, i] <- x[[2L]][rowsL, cols]
        }
        return(rval)
      }
    
      ## compute results
      grid <- sort(unique(c(unlist(lapply(cv, function(x) x[[1]][1, ])))))
      rval <- list(grid = grid)
     
      ## make a matrix with the 'loss' for each fold
      ## at each cp in 'grid'
      
      cv <- lapply(cv, function(x) {
        x[[2]][is.nan(x[[2]]) | is.infinite(x[[2]])] <- NA;
        return(x)
      })
      
      ## oob-loss (computes the average loss)
      oobLoss <- t(sapply(cv, getVals, grid = grid, rowsG = 1L, rowsL = 1L))
      if (length(grid) == 1L) oobLoss <- t(oobLoss)
       
      oobWeights <- matrix(rep(weights, ncol(foldsMat)), ncol = ncol(foldsMat))
      if (attr(foldsMat, "value") == "weights") {
        oobWeights <- oobWeights - foldsMat
      } else {
        oobWeights[foldsMat > 0] <- 0
      }
      oobLoss <- oobLoss / colSums(oobWeights)[!fails]
      rownames(oobLoss) <- paste("fold", 1L:length(cv))
      colnames(oobLoss) <- if (length(grid) > 1L) {
          cn <- c(paste("<=", round(grid[2L], 2)), paste(">", round(grid[-1L], 2)))
      } else {
          cn <- "0"
      }
      
      rval <- append(rval, list(oobloss = oobLoss))
      meanLoss <- colMeans(rval$oobloss, na.rm = TRUE)

      ## cp with minimal loss
      minSubs <- which(meanLoss == min(meanLoss))
      if (length(minSubs) > 1L) minSubs <- max(minSubs)
      rval$cp.hat <- max(0, mean(c(rval$grid, Inf)[minSubs:(minSubs + 1)]))

      rval$foldsMat <- foldsMat
      class(rval) <- "cvloss.tvcm"     
    }

    rval$call <- getCall(object)
    environment(rval$call) <- NULL
    
    if (original) {
      tree$info$cv <- rval
      rval <- tree
    }
    
  } else if (type == "forest") {
    
    rval <- vector(mode = "list", length = 5L)
    names(rval) <- c("node", "coefficients", "contrasts", "error", "folds")
    rval$error$which <-
      which(sapply(cv, function(x) inherits(x, "try-error")))
    rval$error$message <- unlist(cv[rval$error$which])
    cv[rval$error$which] <- NULL
    rval$node <- lapply(cv, function(x) x[[1L]])
    rval$coefficients <- lapply(cv, function(x) x[[2L]])
    rval$contrasts <- lapply(cv, function(x) x[[3L]])   
    rval$folds <- foldsMat
  }
  
  if (control$verbose) cat("\n")
  
  return(rval)
}


print.cvloss.tvcm <- function(x, ...) {
  cat("Cross-Validated Loss", "\n\n")
  cat("Call: ", deparseCall(x$call), "\n\n")
  rval <- colMeans(x$oobloss, na.rm = TRUE)
  if (length(rval) > 10L)
    rval <- rval[seq(1L, length(rval), length.out = 10)]
  print(rval)
  cat("\n")
  cat("cp with minimal oob-loss:", format(x$cp.hat, digits = 3), "\n")
  return(invisible(x))
}


plot.cvloss.tvcm <- function(x, legend = TRUE, details = TRUE, ...) {
  
  xlab <- "cp"
  ylab <- "average loss(cp)"
  type <- "s"
  lpos <- "topleft"
  lsubs <- if (details) 1L:2L else 1L
  if (x$cp.hat < Inf) lsubs <- c(lsubs, 3L)
  lcol <- c("black", "grey80", "black")
  llty <- c(1, 1, 2)

  dotList <- list(...)
  defArgs <- list(type = type, xlab = xlab, ylab = ylab,
                  col = lcol[1], lty = llty[1])
  
  xx <- matrix(x$grid, length(x$grid), nrow(x$oobloss))  
  xx <- rbind(xx, 1.1 * xx[nrow(xx),])

  yy2 <- t(cbind(x$oobloss, x$oobloss[,ncol(x$oobloss)]))
  yy1 <- rowMeans(yy2, na.rm = TRUE)
  
  defArgs$ylim <- range(if (details) yy2 else yy1, na.rm = TRUE)
  
  ## set plot arguments
  call <- list(name = as.name("plot"), x = quote(xx[, 1]), y = quote(yy1))
  call <- appendDefArgs(call, list(...))
  call <- appendDefArgs(call, defArgs)
  type <- call$type
  call$type <- "n"
  llty[1L] <- call$lty
  lcol[1L] <- call$col
  mode(call) <- "call"
  
  ## plot mean validated prediction error
  eval(call)
  
  ## plot details
  if (details)
    matplot(x = xx, yy2, type = type, lty = llty[2L], col = lcol[2L], add = TRUE)

  ## plot average oob-loss as last
  call$name <- as.name("points")
  call$type <- type
  eval(call)
  
  if (legend) {
    ltext <- c("average oob-loss",
               "foldwise oob-loss",
               "cp with minimal oob-loss")
    legend(lpos, ltext[lsubs], col = lcol[lsubs], lty = llty[lsubs])
  }
  
  if (x$cp.hat < Inf) {
    subsMin <- max(which(x$grid <= x$cp.hat))
    minLoss <- colMeans(x$oobloss, na.rm = TRUE)[subsMin]
    segments(x$cp.hat, par()$usr[3], x$cp.hat, minLoss, lty = 2)
    axis(1, x$cp.hat, format(x$cp.hat, digits = 3),
         line = 1, tick = FALSE)
  } 
}
