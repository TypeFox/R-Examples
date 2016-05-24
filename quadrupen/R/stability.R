##' Stability selection for a quadrupen fit.
##'
##' Compute the stability path of a (possibly randomized) fitting
##' procedure as introduced by Meinshausen and Buhlmann (2010).
##'
##' @param x matrix of features, possibly sparsely encoded
##' (experimental). Do NOT include intercept.
##'
##' @param y response vector.
##'
##' @param penalty a string for the fitting procedure used for
##' cross-validation. Either \code{\link{elastic.net}} or
##' \code{"bounded.reg"}.
##'
##' @param subsamples integer indicating the number of subsamplings
##' used to estimate the selection probabilities. Default is 100.
##'
##' @param sample.size integer indicating the size of each subsamples.
##' Default is \code{floor(n/2)}.
##'
##' @param randomize Should a randomized version of the fitting
##' procedure by used? Default is \code{TRUE}. See details below.
##'
##' @param weakness Coefficient used for randomizing. Default is
##' \code{0.5}. Ignored when \code{randomized} is \code{FALSE}. See
##' details below.
##'
##' @param folds list with \code{subsamples} entries with vectors
##' describing the folds to use for the stability procedure. By
##' default, the folds are randomly sampled with the specified
##' \code{subsamples} argument.
##'
##' @param verbose logical; indicates if the progression should be
##' displayed. Default is \code{TRUE}.
##'
##' @param mc.cores the number of cores to use. The default uses all
##' the cores available.
##'
##' @param ... additional parameters to overwrite the defaults of the
##' fitting procedure. See the corresponding documentation
##' (\code{\link{elastic.net}} or \code{\link{bounded.reg}})
##'
##' @return An object of class \code{\linkS4class{stability.path}}.
##'
##' @note When \code{randomized = TRUE}, the \code{penscale} argument
##' that weights the penalty tuned by \eqn{\lambda_1}{lambda1} is
##' perturbed (divided) for each subsample by a random variable
##' uniformly distributed on
##' \if{latex}{\eqn{[\alpha,1]}}\if{html}{[&#945;,1]}\if{text}{\eqn{[alpha,1]}},
##' where
##' \if{latex}{\eqn{\alpha}}\if{html}{&#945;}\if{text}{\eqn{alpha}} is
##' the weakness parameter.
##'
##' If the user runs the fitting method with option
##' \code{'bulletproof'} set to \code{FALSE}, the algorithm may stop
##' at an early stage of the path. Early stops of the underlying
##' fitting function are handled internally, in the following way: we
##' chose to simply skip the results associated with such runs, in
##' order not to bias the stability selection procedure. If it occurs
##' too often, a warning is sent to the user, in which case you should
##' reconsider the grid of \code{lambda1} for stability selection. If
##' \code{bulletproof} is \code{TRUE} (the default), there is nothing
##' to worry about, except a possible slow down when any switching to
##' the proximal algorithm is required.
##'
##' @references N. Meinshausen and P. Buhlmann (2010). Stability
##' Selection, JRSS(B).
##'
##' @seealso \code{\linkS4class{stability.path}} and
##' \code{\link{plot,stability.path-method}}.
##'
##' @examples \dontrun{
##' ## Simulating multivariate Gaussian with blockwise correlation
##' ## and piecewise constant vector of parameters
##' beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
##' Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
##' Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
##' Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
##' diag(Sigma) <- 1
##' n <- 100
##' x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
##' y <- 10 + x %*% beta + rnorm(n,0,10)
##'
##' ## Build a vector of label for true nonzeros
##' labels <- rep("irrelevant", length(beta))
##' labels[beta != 0] <- c("relevant")
##' labels <- factor(labels, ordered=TRUE, levels=c("relevant","irrelevant"))
##'
##' ## Call to stability selection function, 200 subsampling
##' stab <- stability(x,y, subsamples=200, lambda2=1, min.ratio=1e-2)
##' ## Recover the selected variables for a given cutoff
##' ## and per-family error rate, without producing any plot
##' stabpath <- plot(stab, cutoff=0.75, PFER=1, plot=FALSE)
##'
##' cat("\nFalse positives for the randomized Elastic-net with stability selection: ",
##'      sum(labels[stabpath$selected] != "relevant"))
##' cat("\nDONE.\n")
##'}
##' @keywords models, regression
##'
##' @name stability
##' @aliases stability
##' @rdname stability
##'
##' @export
stability <- function(x,
                      y,
                      penalty = c("elastic.net", "bounded.reg"),
                      subsamples  = 100,
                      sample.size = floor(n/2),
                      randomize   = TRUE,
                      weakness    = 0.5,
                      verbose     = TRUE,
                      folds       = replicate(subsamples, sample(1:nrow(x), sample.size), simplify=FALSE),
                      mc.cores    = detectCores(),
                      ...) {

  ## =============================================================
  ## INITIALIZATION & PARAMETERS RECOVERY
  if (Sys.info()[['sysname']] == "Windows") {
    warning("\nWindows does not support fork, enforcing mc.cores to '1'.")
    mc.cores <- 1
  }
  penalty <- match.arg(penalty)
  get.lambda1 <- switch(penalty,
                        "elastic.net" = get.lambda1.l1,
                        "bounded.reg" = get.lambda1.li)
  p <- ncol(x)
  n <- nrow(x)
  user <- list(...)
  defs <- default.args(penalty,nrow(x),ncol(x),user)
  args <- modifyList(defs, user)
  ## overwrite parameters irrelevant in and resampling context
  args$control$verbose <- 0
  args$max.feat <- p
  ## Compute a grid of lambda1 (the smae for each fold)
  if (is.null(args$lambda1)) {
    input <- standardize(x,y,args$intercept,args$normalize,args$penscale)
    args$lambda1 <- get.lambda1(input$xty,args$nlambda1,args$min.ratio)
    rm(input)
  }
  nlambda1 <- length(args$lambda1)
  ## record the basic penscale value for randomizing
  penscale <- args$penscale
  ## Prepare blocs of sub samples to run jobs parallely
  blocs <- split(1:subsamples, 1:mc.cores)

  if (verbose) {
    cat(paste("\n\nSTABILITY SELECTION ",ifelse(randomize,"with","without")," randomization (weakness = ",weakness,")",sep=""))
    cat(paste("\nFitting procedure:",penalty," with lambda2 = ",args$lambda2," and an ",nlambda1,"-dimensional grid of lambda1.", sep=""))
    cat("\nRunning",length(blocs),"jobs parallely (1 per core)")
    cat("\nApprox.", length(blocs[[1]]),"subsamplings for each job for a total of",subsamples)
  }

  ## function to run on each core
  bloc.stability <- function(subsets) {
    select <- Matrix(0,length(args$lambda1),p)
    subsamples.ok <- 0
    for (s in 1:length(subsets)) {
      if (randomize) {args$penscale <- penscale / runif(p,weakness,1)}
      active <- do.call(quadrupen, c(list(x=x[folds[[subsets[s]]], ], y=y[folds[[subsets[s]]]]), args))@active.set
      if (nrow(active) == length(args$lambda1)) {
        subsamples.ok <- subsamples.ok + 1
        select <- select + active
      }
    }
    if (subsamples.ok < 0.5*length(subsets)) {
      cat("\nWarning: more than 50% of the run were discarded in that core due to early stops of the fitting procedure. You should consider largest 'min.ratio' or strongest 'lambda2'.")
    }
    return(select/(subsamples.ok*length(blocs)))
  }

  ## Now launch the B jobs...
  prob.bloc <- mclapply(blocs, bloc.stability, mc.cores=mc.cores)

  ## Construct the probability path
  path <- Matrix(0,nlambda1,p)
  for (b in 1:length(prob.bloc)) {
    path <- path + prob.bloc[[b]]
  }

  return(new("stability.path",
             probabilities = path        ,
             penalty       = penalty     ,
             naive         = args$naive  ,
             lambda1       = args$lambda1,
             lambda2       = args$lambda2,
             folds         = folds       ))

}
