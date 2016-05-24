#' @title Utilities for ``LICORS'' class
#' @name mixed_LICORS-utils
#' @aliases plot.mixed_LICORS summary.mixed_LICORS
#' @description 
#' 
#' The \code{"mixed_LICORS"} class is the objectput from the 
#' \code{\link{mixed_LICORS}} estimator.
#'
NULL

#' @rdname mixed_LICORS-utils
#'
#' @description 
#' \code{plot.mixed_LICORS} gives a visual summary of the estimates
#' such as marginal state probabilities, conditional state probabilities 
#' (= weight matrix), predictive state densities, trace plots for 
#' log-likelihood/loss/penalty.
#'
#' @param x object of class \code{"mixed_LICORS"} 
#' @param type should only \code{"training"}, \code{"test"}, or \code{"both"} be
#' plotted. Default: \code{"both"}.
#' @param cex.axis The magnification to be used for axis annotation relative to 
#' the current setting of \code{cex}.
#' @param cex.lab The magnification to be used for x and y labels 
#' relative to the current setting of \code{cex}.
#' @param cex.main The magnification to be used for main titles relative to the 
#' current setting of \code{cex}.
#' @param line on which margin line should the labels be ploted, 
#' starting at 0 counting objectwards (see also \code{\link[graphics]{mtext}}).
#' @param ... optional arguments passed to \code{plot}, \code{summary},
#' or \code{predict}
#' @keywords hplot
#' @import fields
#' @method plot mixed_LICORS
#' @export
#' @examples
#' # see examples of LICORS-package
#'

plot.mixed_LICORS <- function(x, type = "both", cex.axis = 1.5, cex.lab = 1.5, 
                              cex.main = 2, line = 1.5, ...){
  object <- x
  op <- par(no.readonly = TRUE)
  if (type == "train"){
    state.probs.matrix <- object$conditional.state.probs$opt[object$train.index, ]
    loglik <- object$loglik.trace[, c("train", "train.weighted")]
    mse <- object$MSE.trace[, c("train", "train.weighted")]
  } else if (type == "test"){
    state.probs.matrix <- object$conditional.state.probs$opt[-object$train.index,]
    loglik <- object$loglik.trace[, c("test", "test.weighted")]
    mse <- object$MSE.trace[, c("test", "test.weighted")]
  } else {
    state.probs.matrix <- object$conditional.state.probs$opt
    loglik <- object$loglik.trace
    mse <- object$MSE.trace
  }
  
  state.probs.matrix <- as.matrix(state.probs.matrix)
  states.by.size <- order(colSums(state.probs.matrix), decreasing = TRUE)
  state.probs.matrix <- state.probs.matrix[, states.by.size]
  state.probs <- estimate_state_probs(state.probs.matrix)
  names(state.probs) <- NULL
  
  num.sampless <- nrow(state.probs.matrix)
  
  states.over.time <- rep(NA, object$num.iter['final'])
  states.over.time[1] <- object$num.states["start"]
  
  for (ii in seq_len(object$num.merges)) {
    states.over.time[object$merging.iter[ii]] <- object$num.states["start"] - ii
  }  
  
  FLCs.train.order <- order(rbind(object$LCs[["FLC"]])[object$train.index, ])
  FLCs.train.ordered <- object$LCs[["FLC"]][object$train.index,][FLCs.train.order]

  FLC.pdfs.train <- 
    estimate_LC_pdfs(LCs = cbind(object$LCs[["FLC"]][object$train.index,]), 
                     weight.matrix = object$conditional.state.probs$opt[object$train.index, ][, states.by.size],
                     method = object$control$estimation.method[["FLC"]])
  
  FLC.pdfs.train <- FLC.pdfs.train[FLCs.train.order,]
  
  image.colors <- two.colors(n = 100, "darkblue", "red", "green")
  state.colors <- colorRampPalette(brewer.pal(9, name="Set1"))(object$num.states['opt'])
  
  state.ticks <- pretty(seq_len(object$num.states['opt']))
  state.probs.ticks <- pretty(state.probs)
  LC.ticks <- pretty(seq_len(num.sampless), 5)
  iter.ticks <- pretty(seq_len(object$num.iter['final']), 
                           object$num.iter['final'])
  FLC.ticks <- pretty(object$LCs[["FLC"]][object$train.index,])
  pdf.FLC.ticks <- pretty(FLC.pdfs.train)
  ###########################
  ### Start the plot
  ###########################
  
  par(mfrow = c(1, 1), mar = c(5, 5, 2, 5), 
      cex.axis = cex.axis, cex.lab = cex.lab)
  layout(matrix(c(1, 2, 3, 3, 4, 5, 6, 7), byrow = FALSE, ncol = 2))
  # conditional FLC distributions
  par(mar = c(1, 4, 1, 1))
  matplot(FLCs.train.ordered, FLC.pdfs.train, type= "l", lwd = 2, lty = 1,
          col = state.colors, axes = FALSE, ylab = "", xlab = "", 
          main = "", cex.main = cex.main)
  box() 
  axis(1, at = FLC.ticks, labels = paste(FLC.ticks))
  axis(4, at = pdf.FLC.ticks, labels = paste(pdf.FLC.ticks))
  
  abline(v = FLCs.train.ordered[apply(FLC.pdfs.train, 2, which.max)], 
         col = state.colors, lty =2)
  mtext("p(x|S)", 2, line = line, cex = cex.lab)
  legend("topright", ifelse(type == "both", "training & test", "training only"))
  
  # distribution over states
  par(mar = c(1,4,1,1))
  barplot(state.probs, space = 0, main = "" , xlab = "", ylab = "", 
          col = state.colors, 
          cex.axis = cex.axis, cex.main = cex.main, 
          ylim = c(0, max(state.probs) * 1.05),
          axes = FALSE, xlim = c(0 + 0.25, object$num.states['opt'] - 0.25))
  legend("topright", ifelse(type == "both", "training & test", type))
  
  mtext("P(S)", 2, line = line, cex = cex.lab)
  box()
  axis(4, at = state.probs.ticks, labels = paste(state.probs.ticks))
  
  par(mar = c(5, 4, 0.1, 1))
  image2(state.probs.matrix, axes = FALSE, legend = FALSE,
        col = image.colors, xlab = "", ylab = "", main = "", cex.main = cex.main)
  abline(v = 0:object$num.states['opt'] + 0.5, col = "gray", lwd = 1)
  if (type == "both"){
    abline(h = length(object$theta) - length(-object$train.index), col = "white", lwd = 2)
  }
  box(col = "black")
  mtext("state id", 1, line = line + 2, cex = cex.lab)
  mtext("P(S|x, PLC)", 2, line = line, cex = cex.lab)
  mtext("PLC id", 4, line = line + 2, cex = cex.lab, 
        at = nrow(state.probs.matrix) / 3)
  
  if (type == "both"){
    mtext("test", 2, line = 0.25, at = (length(object$theta) - length(-object$train.index))*0.5)
    mtext("training", 2, line = 0.25, at = nrow(state.probs.matrix)*4/5 )
  }
  
  axis(1, at = state.ticks, labels = paste(state.ticks), 
       col = state.colors)
  axis(4, at = LC.ticks, labels = rev(paste(LC.ticks)))
  
  
  # Right column of the plot
  par(mar = c(0, 5, 1, 3))
  if (type != "both"){
    ts.plot(ts(loglik), pch = 1, main = "", ylab = "", xlab = "", axes = FALSE)
    abline(v = c(0, object$merging.iter - 0.5))
    axis(2)
    axis(1, at = iter.ticks, labels = paste(iter.ticks))
    box()
    mtext("Iteration", 1, line = line, cex = cex.lab)
    mtext(paste0("Log-likelihood (", type, ")"), 2, line = line, cex = cex.lab)
    
    par(new = TRUE)
    plot(states.over.time, pch = 19, xlab = "", ylab = "", axes = FALSE)
    states.over.time = na.locf(states.over.time)
    lines(states.over.time, lwd=2)
    axis(2, at = pretty(object$num.states["start"]:object$num.states["final"], 
                        n = object$num.merges))
    mtext("Number of states", 4, line = 2, cex = cex.lab)
    
  } else {
    matplot(loglik, ylab = "", main = "", type = "l", lty = c(1, 2, 1, 2), 
            col = c(1, 1, 2, 2), lwd = 2, axes = FALSE, cex.main = cex.main)
    box()
    axis(4)
    abline(v = object$num.iter['opt'], lty = 3, lwd = 2)
    mtext("log-likelihood", 2, line = line, cex = cex.lab)
  }
  
  par(mar = c(0,5,0,3))
  # MSE comparison
  matplot(mse, ylab = "", main = "", type = "l", lty = c(1, 2, 1, 2), 
          col = c(1, 1, 2, 2), lwd = 2, axes = FALSE, cex.main = cex.main)
  box()  
  axis(4)
  
  abline(v = object$num.iter['opt'], lty = 3, lwd = 2)
  mtext("MSE", 2, line = line, cex = cex.lab)
  
  par(mar = c(4,5,0,3))
  matplot(object$penalty[1:object$num.iter['final'], "entropy"], type = "l", 
          lty = 1:2, col = 1:2, lwd = 2, main = "", ylab = "", 
          cex.main = cex.main, axes = FALSE)
  box()
  axis(1, at = iter.ticks, labels = paste(iter.ticks))
  axis(4)

  abline(v = object$num.iter['opt'], lty = 3, lwd = 2)
  mtext("Iteration", 1, line = line + 1, cex = cex.lab)  
  mtext("penalty", 2, line = line, cex = cex.lab)
  
  par(mar = c(1,7,1,3))
  plot.new()
  legend("topright", colnames(mse), lty = c(1, 2, 1, 2), col = c(1, 1, 2, 2), 
         lwd=2, cex = cex.lab)
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  rect(0, seq(0, 1, length = 100)[-100],
       .1, seq(0, 1, length = 100)[-1],
       col = image.colors, border = NA)

  axis(2, at = pretty(seq(0, 1, length = 100)), las = 0)
  mtext("probability", side = 2, cex = cex.lab*0.75, line = 2, at = 0.25)
  
  par(op)
}

#' @rdname mixed_LICORS-utils
#' @description 
#' \code{summary.mixed_LICORS} prints object a summary of the estimated LICORS model.
#'
#' @param object object of class \code{"mixed_LICORS"} 
#' @keywords model nonparametric
#' @export
#' @method summary mixed_LICORS
#' @examples
#' # see examples in LICORS-package
#'

summary.mixed_LICORS <- function(object, ...){
  
  cat(rep("*", 20))
  cat("\n")
  
  if (object$converged){
    cat("Mixed LICORS converged after", object$num.iter["final"], "iterations
        (", object$num.states["final"], "states).")
  } else {
    cat("Mixed LICORS did NOT converge. \n")
    cat("Check trace plots to see if the solution is an interior optimum (vertical dashed bar in trace plots).\n \n")
  }
  cat("\n")
  cat("The best model W* was achieved at iteration", object$num.iter['opt'], 
      "with", object$num.states["start"], 
      "remaining states (starting from", object$num.states["start"],"). \n")
    
  entropies <- rep(NA, 3)
  names(entropies) <- c("both", "training", "test")
  entropies["both"] <- compute_mixture_penalty(object$conditional.state.probs$opt, 
                                               type = "entropy", base = "num.states")
  entropies["training"] <- 
    compute_mixture_penalty(object$conditional.state.probs$opt[object$train.index,], 
                            type = "entropy", base = "num.states")
  entropies["test"] <- 
    compute_mixture_penalty(object$conditional.state.probs$opt[-object$train.index,], 
                             type = "entropy", base = "num.states")
  cat("\n")  
  cat("The average entropy penalty for W* equals", 
      round(object$penalty[object$num.iter['opt'],"entropy"]*100, 1), 
      "%, where 0% is unique cluster assignment, and 100% is uniformly at random.\n\n") 
  cat("Entropy for both vs training vs test (in %): \n")
  cat(round(cbind(entropies * 100), 1))
  cat("\n")
  cat(rep("*", 20))
}

#' @rdname mixed_LICORS-utils
#' @description 
#' \code{predict.mixed_LICORS} predicts FLCs based on PLCs given a fitted
#' mixed LICORS model.
#' This can be done on an iterative basis, or for a selection of future PLCs.
#'
#' @param new.LCs a list with PLC configurations to predict FLCs given these PLCs
#' @keywords model nonparametric
#' @export
#' @method predict mixed_LICORS
#' @examples
#' # see examples in LICORS-package
#'

predict.mixed_LICORS <- function(object, new.LCs = list(PLC = NULL), ...) {
  
  prediction <- list(FLC = NULL,
                     states = NULL,
                     weight.matrix = NULL)
  
  if (is.null(new.LCs$PLC)) {
    # predict the model fit on the historical data
    prediction$weight.matrix <- object$conditional.state.probs$opt
  } else {
    pdfs.of.PLC.new <- estimate_LC_pdfs(LCs = object$LCs$PLC, 
                                        weight.matrix = 
                                          object$conditional.state.probs$opt,
                                        method = 
                                          object$control$estimation.method[["PLC"]],
                                        eval.LCs = new.LCs$PLC)
    prediction$weight.matrix <- 
      estimate_state_probs(weight.matrix = object$conditional.state.probs$opt,
                           states = NULL, 
                           pdfs = list(PLC = pdfs.of.PLC.new, 
                                       FLC = NULL),
                           num.states = ncol(pdfs.of.PLC.new))
  }
  
  prediction$state <- weight_matrix2states(prediction$weight.matrix)
  prediction$FLC <- 
    predict_FLC_given_PLC(train = list(data = list(FLC = object$LCs$FLC,
                                                   PLC = object$LCs$FLC),
                                       weight.matrix = 
                                         object$conditional.state.probs$opt),
                          test = list(weight.matrix = prediction$weight.matrix),
                          ...)
  if (is.null(new.LCs$PLC)) {
    prediction$residuals <- object$LCs$FLC - prediction$FLC
  }
  return(prediction)
}

#' @rdname mixed_LICORS-utils
#' @description 
#' \code{complete_LICORS_control} completes the controls for 
#' the mixed LICORS estimator. Entries of the list are:
#' 
#' 'loss' an R function specifying the loss for cross-validation (CV). 
#' Default: mean squared error (MSE), i.e. \code{loss = function(x, xhat) mean((x-xhat)^2)}
#'
#' 'method' a list of length \eqn{2} with arguments \code{PLC} and 
#' \code{FLC} for the method of density estimation in each 
#' (either \code{"normal"} or \code{"nonparametric"}).
#' 
#' 'max.iter' maximum number of iterations in the EM
#' 
#' 'trace' if > 0 it prints output in the console as the EM is running
#' 
#' 'sparsity' what type of sparsity (currently not implemented)
#' 
#' 'lambda' penalization parameter; larger lambda gives sparser weights
#' 
#' 'alpha' significance level to stop testing. Default: \code{alpha = 0.01}
#' 
#' 'seed' set seed for reproducibility. Default: \code{NULL}. If 
#' \code{NULL} it sets a random seed and then returns this seed in the output.
#' 
#' 'CV.train.ratio' how much of the data should be training data.
#' Default: \code{0.75}, i.e., \eqn{75\%} of data is for training
#' 
#' 'CV.split.random' logical; if \code{TRUE} training and test data
#' are split randomly; if \code{FALSE} (default) it uses the first
#' part (in time) as training, rest as test.
#' 
#' 'estimation' a list of length \eqn{2} with arguments \code{PLC} and 
#' \code{FLC} for the method of density estimation in each 
#' (either \code{"normal"} or \code{"nonparametric"}).
#' 
#' @param control a list of controls for \code{"mixed_LICORS"}.
#' @keywords model nonparametric
#' @export
#' @examples
#' # see examples in LICORS-package
#'

complete_LICORS_control <- 
  function(control = list(alpha = 0.01, 
                          CV.split.random = FALSE,
                          CV.train.ratio = 0.75,
                          lambda = 0,
                          max.iter = 500, 
                          seed = NULL, 
                          sparsity = "stochastic",
                          trace = 0,
                          loss = function(x, xhat) 
                            mean((x - xhat)^2),
                          estimation.method = 
                            list(PLC = "normal", FLC = "nonparametric"))) {

  # complete control settings
  if (is.null(control$loss)) {
    control$loss <- function(x, xhat) mean((x - xhat)^2)
  }
  if (is.null(control$estimation.method)) {
    control$estimation.method <- 
      list(PLC = "normal", FLC = "nonparametric")
  }
  if (is.null(control$alpha)) {
    control$alpha <- 0.01
  }
  if (is.null(control$CV.split.random)) {
    control$CV.split.random <- FALSE
  }
  if (is.null(control$CV.train.ratio)) {
    control$CV.train.ratio <- 0.75
  }
  if (is.null(control$lambda)) {
    control$lambda <- 0
  }
  if (is.null(control$max.iter)) {
    control$max.iter <- 500
  }
  if (is.null(control$sparsity)) {
    control$sparsity <- "stochastic"
  }
  if (is.null(control$trace)) {
    control$trace <- 0
  }
  if (is.null(control$seed)){
    control$seed <- sample.int(10^6, 1)
  }
  return(control)
}

