#' @title Estimate PLC/FLC distributions for all states
#' @aliases estimate_LC.pdf.state
#' @description  
#' \code{\link{estimate_LC_pdfs}} estimates the PLC and FLC distributions for 
#' each state \eqn{k = 1, \ldots, K}. It iteratively applies 
#' \code{\link{estimate_LC.pdf.state}}.
#' 
#' \code{\link{estimate_LC.pdf.state}} estimates the PLC and FLC 
#' distributions using weighted maximum likelihood (\code{\link[stats]{cov.wt}}) 
#' and nonparametric kernel density estimation (\code{\link{wKDE}}) for one
#'  (!) state.
#' 
#' @param LCs matrix of PLCs/FLCs. This matrix has \eqn{N} rows and 
#' \eqn{n_p} or \eqn{n_f} columns (depending on the PLC/FLC dimensionality)
#' @param weight.matrix \eqn{N \times K} weight matrix
#' @param states vector of length \eqn{N} with entry \eqn{i} being the label 
#' \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param method type of estimation: either a (multivariate) Normal distribution 
#' (\code{"normal"}) or nonparametric with a 
#' kernel density estimator (\code{method = "nonparametric"}).  
#' For multivariate distributions (as usual for PLCs) only
#' \code{'normal'} should be used due to computational efficiency and 
#' statistical accuracy.
#' @param eval.LCs on what LCs should the estimate be evaluated? If \code{NULL} then densities
#' will be evaluated on the training data \code{LCs}
#' @return
#' \code{\link{estimate_LC_pdfs}} returns an \eqn{N \times K} matrix.
#' 
#' @keywords nonparametric multivariate distribution
#' @export
#' @examples
#' set.seed(10)
#' WW = matrix(runif(10000), ncol = 10)
#' WW = normalize(WW)
#' temp_flcs = cbind(sort(rnorm(nrow(WW))))
#' temp_flc_pdfs = estimate_LC_pdfs(temp_flcs, WW)
#' matplot(temp_flcs, temp_flc_pdfs, col = 1:ncol(WW), type = "l", 
#'         xlab = "FLCs", ylab = "pdf", lty = 1)
#' 
#' 
estimate_LC_pdfs <- function(LCs, weight.matrix = NULL,
                             method = 
                               c("nonparametric", "normal", "huge"),
                             eval.LCs = NULL) {
  
  method <- match.arg(method)
  if (is.null(weight.matrix)) {
    weight.matrix <- matrix(rep(1, nrow(LCs)), ncol = 1)
  }
  num.states <- ncol(weight.matrix)
  states <- weight_matrix2states(weight.matrix)
   
  if (is.null(eval.LCs)) {
    LCs <- cbind(LCs)
    num.evals <- nrow(LCs)
  } else {
    eval.LCs <- cbind(eval.LCs)
    num.evals <- nrow(eval.LCs)
  }
  if (is.null(num.evals)) {
    num.evals <- 1
  }
  pdf.all.states <- Matrix(0, ncol = num.states, nrow = num.evals,
                           sparse = TRUE)
  for (ii in seq_len(num.states)) {
    if (all(weight.matrix[, ii] == 0)) {
      # set probability to 0 if all weights are 0
      pdf.all.states[, ii] <- 0
    } else {
      pdf.all.states[, ii] <- 
        estimate_LC_pdf_state(state = ii, 
                              states = states,
                              weights = weight.matrix[, ii],
                              eval.LCs = eval.LCs,
                              LCs = LCs, 
                              method = method)
    }
  }
  invisible(pdf.all.states)
} 


#' @rdname estimate_LC_pdfs
#' @param state integer; which state-conditional density should be 
#' estimated
#' @param weights weights of the samples. Either a i) length \eqn{N} vector with the 
#' weights for each observation; ii) \eqn{N \times K} matrix, where the
#' \code{state} column of that matrix is used as a weight-vector.
#' @keywords nonparametric multivariate distribution
#' @return
#' \code{\link{estimate_LC.pdf.state}} returns a vector of length \eqn{N} 
#' with the state-conditional density evaluated at \code{eval.LCs}.
#' @export
#' @examples
#' 
#' ######################
#' ### one state only ###
#' ######################
#' temp_flcs <- temp_flcs[order(temp_flcs)]
#' temp_flc_pdf = estimate_LC_pdf_state(state = 3, 
#'                                      LCs = temp_flcs, 
#'                                      weights = WW)
#'     
#' plot(temp_flcs, temp_flc_pdf, type = "l", xlab = "FLC", ylab = "pdf")

estimate_LC_pdf_state <- function(state, states = NULL, 
                                  weights = NULL, 
                                  LCs = NULL, 
                                  eval.LCs = NULL, 
                                  method = 
                                    c("nonparametric", "normal", "huge")) {
  
  method <- match.arg(method)
  
  if (is.null(LCs)){
    stop("You have to provide data (LCs) to estimate the densities.")
  } else {
    LCs <- cbind(LCs)
  }
  if (is.null(eval.LCs)) {
    eval.LCs <- LCs
  } else {
    eval.LCs <- cbind(eval.LCs)
  }
  num.states <- ncol(eval.LCs)
  
  if (is.null(states) && is.null(weights)) {
    stop("Either provide state assignment or weights.")
  }
  
  if (is.matrix(weights)){
    # necessary below for optimal bandwidth selection
    states <- weight_matrix2states(weights)
    weights <- weights[, state]
    if (any(is.na(weights / sum(weights)))) {
      stop(paste0("All weights are 0; normalizing introduces 'NA'. 
                   Please remove state", state, "."))
    }
  } 
  
  in.current.state <- which(states == state)
  LCs.in.current.state <- LCs[in.current.state, ]
  if (is.null(weights)) {
    effective.sample.size <- length(in.current.state)
  } else {
    effective.sample.size <- sum(weights)
  }
  
  # state adaptive bandwidth selection adjust the bw.nrd0() rule to only 
  # observations in the given state
  if (length(in.current.state) < 5){
    in.current.state <- which(weights >= quantile(weights, 0.95))
  }
  
  if (length(in.current.state) < 5) {
    optimal.bw <- 0.1
  } else {
    optimal.bw <- bw.nrd0(c(LCs)[in.current.state])
  }
  
  if (num.states == 1) {
    # univariate density estimation
    switch(method,
           nonparametric = {
             if (is.null(weights) || is.na(weights)) {
               LC.pdf.state <- wKDE(LCs.in.current.state, 
                                    eval.points = eval.LCs, 
                                    weights = NULL,
                                    bw = optimal.bw)
             } else {
               LC.pdf.state <- wKDE(LCs, eval.points = eval.LCs, 
                                    weights = weights / sum(weights),
                                    bw = optimal.bw)
             }
           },
           normal = {
             if (is.null(weights) || is.na(weights)) {
               mu <- mean(LCs.in.current.state)
               sigma <- sd(c(LCs.in.current.state))
             } else {
               temp <- cov.wt(LCs, wt = weights, method = "ML")
               mu <- temp$center
               sigma <- sqrt(temp$cov)
               rm(temp)
             }
             if (effective.sample.size < 5) {
               sigma <- apply(LCs, 2, sd)
             }
             LC.pdf.state <- dnorm(c(eval.LCs), mean = mu, sd = sigma)
           })
  } else {  # multivariate density estimation
    switch(method,
           nonparametric = {
             stop('Multivariate nonparametrc not implemented yet.')
             if (is.null(weights) || is.na(weights)) {
               LC.pdf.state <- mv_wKDE(LCs.in.current.state, 
                                       eval.points = eval.LCs)
             } else {
               LC.pdf.state <- mv_wKDE(LCs, eval.points = eval.LCs, 
                                       weights = weights / sum(weights))
             }
           },
           normal = {
             if (effective.sample.size < 5) {
               Sigma.mat <- cov(LCs)
               mu.vec <- apply(LCs, 2, median)
             } else {
               if (is.null(weights) || is.na(weights)) {
                 mu.vec <- colMeans(LCs.in.current.state)
                 Sigma.mat <- cov(LCs.in.current.state)
               } else {
                 temp <- cov.wt(LCs, wt = weights, method = "ML")
                 mu.vec <- temp$center
                 Sigma.mat <- temp$cov
                 rm(temp)
               }
             }
             Sigma.mat <- Sigma.mat + 
               1/100 * mean(diag(Sigma.mat)) * diag(1, ncol(Sigma.mat))
             LC.pdf.state <- dmvnorm(eval.LCs, mean = mu.vec, sigma = Sigma.mat)
           },
           huge = {
             glasso.est <- huge.select(huge(LCs.in.current.state, 
                                            cov.output = TRUE, method = "glasso", 
                                            verbose = FALSE), 
                                       verbose = FALSE)
             LC.pdf.state <- dmvnorm(eval.LCs, 
                       mean = apply(LCs.in.current.state, 2, median), 
                       sigma = as.matrix(glasso.est$opt.cov))
           })
  }
  LC.pdf.state[is.na(LC.pdf.state)] <- .Machine$double.eps^0.5
  invisible(LC.pdf.state)
}
