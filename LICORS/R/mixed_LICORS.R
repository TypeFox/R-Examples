#' @title Mixed LICORS: An EM-like Algorithm for Predictive State Space Estimation
#'
#' @description \code{mixed_LICORS} is the core function of this
#' package as it estimates the ``parameters'' in the model for the
#' spatio-temporal process.  \deqn{ P(X_1, \ldots, X_{\tilde{N}})
#' \propto \prod_{i=1}^{N} P(X_i \mid \ell^{-}_i) = \prod_{i=1}^{N}
#' P(X_i \mid \epsilon(\ell^{-}_i)) .  }
#' 
#' @param LCs list of PLCs and FLCs matrices (see output of
#' \code{\link{data2LCs}} for details and formatting).
#' @param num.states.init number of states to start the EM algorithm
#' @param initialization a a) character string, b) vector, or c)
#' matrix. a) results \code{num.states.init} many states initialized
#' by passing the character string as \code{method} argument of
#' \code{\link{initialize_states}}; if b) the vector will be taken as
#' initial state labels; if c) the matrix will be taken as initial
#' weights. Note that for both b) and c) \code{num.states.init} will
#' be ignored.  \eqn{k = 1, \ldots, K} of PLC \eqn{i}
#' @param control a list of control settings for the EM algorithm. 
#' See \code{\link{complete_LICORS_control}} for details.
#' @return
#' An object of class \code{"LICORS"}.
#' @keywords nonparametric cluster multivariate distribution
#' @export
#' @seealso \code{\link{plot.mixed_LICORS}},
#' \code{\link{summary.mixed_LICORS}}
#' @examples
#' \dontrun{
#' data(contCA00)
#'  
#' LC_geom = setup_LC_geometry(speed = 1,
#'                             horizon = list(PLC = 2, FLC = 0),
#'                             shape ="cone")
#' bb = data2LCs(t(contCA00$observed), LC.coordinates = LC_geom$coordinates)
#' 
#' mm = mixed_LICORS(bb, num.states.init = 15, init = "KmeansPLC", 
#'                   control = list(max.iter = 50, lambda = 0.001))
#' plot(mm)
#' ff_new = estimate_LC_pdfs(bb$FLC,
#'                           weight.matrix = mm$conditional_state_probs, 
#'                           method = "nonparametric")
#' matplot(bb$FLC, ff_new, pch = ".", cex = 2)
#'}
#'

mixed_LICORS <- function(LCs = list(PLC = NULL, FLC = NULL,
                                    dim = list(original = NULL, 
                                               truncated = NULL)),
                         num.states.init = NULL,
                         initialization = NULL,
                         control = list(max.iter = 500, alpha = 0.01,
                                        trace = 0, lambda = 0,
                                        sparsity = "stochastic",
                                        CV.split.random = FALSE,
                                        CV.train.ratio = 0.75, seed = NULL,
                                        loss = function(x, xhat) 
                                          mean((x - xhat)^2),
                                        estimation.method = 
                                          list(PLC = "normal",
                                               FLC = "nonparametric"))) {
  
  control <- complete_LICORS_control(control)
  set.seed(control$seed)
  if (length(dim(LCs$FLC)) == 0) {
    LCs$FLC <- cbind(LCs$FLC)
  }
  NN <- nrow(LCs$PLC)
  
  if (is.null(LCs$dim)){
    LCs$dim <- list(original = NULL, truncated = c(1, NN))  
  }
  
  if (is.null(initialization) & is.null(num.states.init)){
    stop("You must provide either 
         - an integer specifying the number of states
         - a vector with the initial state space labels
         - or a matrix with the initial weights
         ")
  }
  states <- rep(NA, NN)
  
  if (control$CV.split.random){
    train.index <- sample.int(NN, size = floor(control$CV.train.ratio * NN), 
                              replace = FALSE)
  } else {
    train.index <- seq_len(floor(control$CV.train.ratio * NN))
  }
  
  train.data <- list(PLC = LCs$PLC[train.index, ],
                     FLC = cbind(LCs$FLC[train.index, ]))
  test.data <- list(PLC = LCs$PLC[-train.index, ],
                    FLC = cbind(LCs$FLC[-train.index, ]))
  num.samples <- c("train" = length(train.index),
                   "test" = NN - length(train.index),
                   "all" = NN)
  #####################################################
  #### initialize variables to keep track of EM updates
  #####################################################
  # mixture weight penalty
  penalty <- matrix(NA, nrow = control$max.iter, ncol = 4)
  colnames(penalty) <- c("L2", "entropy", "nec", "nec.weighted")
  
  loglik <- matrix(-Inf, ncol = 5, nrow = control$max.iter)
  colnames(loglik) <- c("train", "train.weighted", "train.weighted.update",
                        "test", "test.weighted")
  MSE <- matrix(NA, ncol = 4, nrow = control$max.iter)
  colnames(MSE) <- c("train", "train.weighted", "test", "test.weighted")

  weight.matrices <- list(train = list(tmp = NULL, final = NULL, opt = NULL),
                          test = list(tmp = NULL, final = NULL, opt = NULL))
  ##########################
  #### Initialize states ###
  ##########################
  if (is.character(initialization)) {   
    if (is.null(num.states.init)){
      stop("You must provide the number of maximum states to start with.")
    } else {
      # assign states by kmeans or randomly
      states[train.index] <-
          initialize_states(num.states = num.states.init, 
                            method = initialization,
                            LCs = train.data)                       
      # state vector initialization for test data is irrelevant; but need
      # the full vector later
      states[-train.index] <-
          initialize_states(num.states = num.states.init,
                            num.samples = num.samples["test"],
                            method = "random")
    }
  } else {
    if (is.vector(initialization)){
      states <- initialization
    } else if (is.matrix(initialization)){
      weight.matrices$train$tmp <- initialization[train.index, ]
      states <- weight_matrix2states(weight.matrices$train$tmp)
    } else {
      stop("Initial states must be either a 
           i) character string describing the initialization method,
           ii) a vector with the labels (length = N)
           iii) or a weight matrix (dim = N x K)")
    }
    
  }
  state.vector <- list(train = states[train.index],
                       test = states[-train.index])
  
  if (is.null(weight.matrices$train$tmp)){
    weight.matrices$train$tmp <-
        states2weight_matrix(states = state.vector$train, 
                             num.states.total = num.states.init)
  }
  num.states <- c(init = max(states), current = max(states),
                  opt = NA, final = NA)

  pdfs <- list()
  pdfs$FLC$train$one.cluster <- 
    estimate_LC_pdf_state(state = 1,
                          states = rep(1, num.samples["train"]),
                          LCs = train.data$FLC,
                          method = control$estimation.method[["FLC"]])
  loglik.one.state <-
      compute_LICORS_loglik(weight.matrix = 1,
                            pdfs.FLC = cbind(pdfs$FLC$train$one.cluster))
    
  penalty[1, "L2"] <- compute_mixture_penalty(weight.matrices$train$tmp,
                                              "Lq", q = 2)
  penalty[1, "entropy"] <- compute_mixture_penalty(weight.matrices$train$tmp,
                                                   "entropy",
                                                   base = "num.states")
  num.states["init"] <- ncol(weight.matrices$train$tmp)

  pdfs$FLC$train <- estimate_LC_pdfs(LCs = train.data$FLC, 
                                     weight.matrix =
                                       weight.matrices$train$tmp,
                                     method = control$estimation.method[["FLC"]])
  pdfs$PLC$train <- estimate_LC_pdfs(LCs = train.data$PLC, 
                                     weight.matrix =
                                       weight.matrices$train$tmp,
                                     method = control$estimation.method[["PLC"]])
  # predictions are stored in a list, not a matrix, since the could be
  # multidimensional themselves
  prediction <- list()
  prediction$FLC$train <-
    predict_FLC_given_PLC(train = list(data = list(FLC = train.data$FLC, 
                                                   PLC = NULL),
                                       weight.matrix = weight.matrices$train$tmp),
                          test = list(weight.matrix = weight.matrices$train$tmp),
                          type = "mean")
  prediction$FLC$train.weighted <-
    predict_FLC_given_PLC(train = list(data = list(FLC = train.data$FLC),
                                       weight.matrix = weight.matrices$train$tmp),
                          test = list(weight.matrix = weight.matrices$train$tmp),
                          type = "weighted.mean")
  
  MSE[1, "train"] <- control$loss(train.data$FLC, prediction$FLC$train)
  MSE[1, "train.weighted"] <- control$loss(train.data$FLC,
                                   prediction$FLC$train.weighted)
  
  pdfs$PLC$test <- estimate_LC_pdfs(LCs = train.data$PLC, 
                                    weight.matrix =
                                      weight.matrices$train$tmp,
                                    method = control$estimation.method[["PLC"]],
                                    eval.LCs = test.data$PLC)
  weight.matrices$test$tmp <-
      estimate_state_probs(weight.matrix = weight.matrices$train$tmp,
                           pdfs = list(PLC = pdfs[["PLC"]][["test"]], FLC = NULL),
                           num.states = num.states["init"])
  
  state.vector$test <- weight_matrix2states(weight.matrices$test$tmp)
  prediction$FLC$test.weighted <-
    predict_FLC_given_PLC(train = list(data = list(FLC = train.data$FLC),
                                       weight.matrix = weight.matrices$train$tmp),
                          test = list(weight.matrix = weight.matrices$test$tmp),
                          type = "weighted.mean")
  prediction$FLC$test <-
    predict_FLC_given_PLC(train = list(data = list(FLC = train.data$FLC),
                                       weight.matrix = weight.matrices$train$tmp),
                          test = list(weight.matrix = weight.matrices$test$tmp),
                          type = "mean")
  
  MSE[1, "test"] <- control$loss(test.data$FLC, prediction$FLC$test)
  MSE[1, "test.weighted"] <- control$loss(test.data$FLC, prediction$FLC$test.weighted)
    
  loglik[1, c("train.weighted", "train", "train.weighted.update")] <- 
    compute_LICORS_loglik(pdfs$FLC$train, 
                          weight.matrices$train$tmp, 
                          lambda = control$lambda)

  pdfs$FLC$test <- estimate_LC_pdfs(LCs = train.data$FLC, 
                                    weight.matrix =
                                      weight.matrices$train$tmp,
                                    method = control$estimation.method[["FLC"]],
                                    eval.LCs = test.data$FLC)
  
  loglik[1, "test.weighted"] <-
      compute_LICORS_loglik(weight.matrices$test$tmp, 
                            pdfs$FLC$test, 
                            lambda = control$lambda)
  loglik[1, "test"] <-
      compute_LICORS_loglik(states2weight_matrix(state.vector$test,
                                                 num.states.total =
                                                   num.states["init"]),
                            pdfs$FLC$test,
                            lambda = control$lambda)
  weight.matrices$train$final <- weight.matrices$train$tmp
  prediction$final <- prediction$FLC$test
  
  num.iter <- c(opt = 1, final = NA, opt = NA, max = control$max.iter)
  
  weight.matrices$train$opt <- weight.matrices$train$tmp
  prediction$opt <- prediction$FLC$test
  prediction$opt.weighted <- prediction$FLC$test
  
  weight.matrices$test$opt <- weight.matrices$test$tmp
 
  stop.iterations <- FALSE
  converged.tmp <- FALSE
  overall.converged <- FALSE
  
  merging.iter = c()
  if (control$trace > 0) {
    cat("Starting EM iterations \n")
  }
  for (ii in 2:control$max.iter){
    if (control$trace > 0) {
      cat("Start iteration", ii, "\n")
    }
    merged <- FALSE
    do.merge.tmp <- FALSE
    more.sparse.tmp <- TRUE
    
    iter.final <- ii
    num.states.tmp <- ncol(weight.matrices$train$tmp)
    # E step
    weight.matrices$train$tmp <- 
      estimate_state_probs(weight.matrix = weight.matrices$train$tmp,
                           pdfs = list(PLC = pdfs$PLC$train, 
                                       FLC = pdfs$FLC$train))
    if (any("deterministic" == control$sparsity) && control$lambda != 0){
      if (!merged && ii > 1000){#} && ii < control$max.iter){
        if (control$trace > 0){
          cat("sparsity enforced! \n")
          cat("Penalty before:",
              round(compute_mixture_penalty(weight.matrices$train$tmp,
                                            "entropy",
                                            base = "num.states") * 100, 1),
              "% \n")
        }
        weight.matrices$train$tmp <-
            sparsify_weights(weight.matrices$train$tmp,
                             NULL, lambda = control$lambda)
  	    if (control$trace > 0){
          cat("Penalty after:",
              round(compute_mixture_penalty(weight.matrices$train$tmp,
                                            "entropy",
                                            base = "num.states") * 100, 1),
              "% \n")
        }
      }
    }

    if (merged) {
      pdfs$PLC$train <- estimate_LC_pdfs(LCs = train.data$PLC, 
                                         weight.matrix =
                                           weight.matrices$train$tmp,
                                         method = control$estimation.method[["PLC"]])
      
      pdfs$FLC$train <- estimate_LC_pdfs(LCs = train.data$FLC, 
                                         weight.matrix =
                                           weight.matrices$train$tmp,
                                         method = control$estimation.method[["FLC"]])
      # M step
      weight.matrices$train$tmp <-
        estimate_state_probs(weight.matrix = weight.matrices$train$tmp,
                             pdfs = list(PLC = pdfs$PLC$train,
                                         FLC = pdfs$FLC$train))
      max.weight.diff <- 1
    } else {  # not merged
      # maximum weight difference
      max.weight.diff <-
          max(sqrt(rowMeans((weight.matrices$train$tmp - weight.matrices$train$final)^2)))
    }

    if (max.weight.diff < 10^(-3)){
      converged.tmp <- TRUE
    }
    
    if (converged.tmp) {
      do.merge.tmp <- TRUE
      if (any("deterministic" == control$sparsity) && control$lambda != 0){
        do.merge.tmp <- FALSE
        penalty.old.tmp <- penalty[ii - 1, "entropy"]
        if (control$trace > 0){
          cat("converged and sparsity enforced! \n")
          cat("Penalty before:", round(penalty.old.tmp * 100, 4), "% \n")
        }
        if (control$lambda > 0) {
          # sparsity
          weight.matrices$train$tmp <- 
            sparsify_weights(weight.matrices$train$tmp, NULL,
                             lambda = control$lambda)        
        }
        penalty.new.tmp <- compute_mixture_penalty(weight.matrices$train$tmp,
                                                    "entropy")

        more.sparse.tmp <- (penalty.old.tmp >= penalty.new.tmp)
        if (control$trace > 0){
          cat("Penalty after:", round(penalty.new.tmp * 100, 4), "% \n")
        }
        
        max.weight.diff <-
          max(sqrt(rowMeans((weight.matrices$train$tmp - weight.matrices$train$final)^2)))
        
        if (max.weight.diff < 10^(-3.01) 
            || ( any( abs(loglik[ii - 1, "train.weighted"] -
                          loglik[ii - seq_len(min(ii - 1, 10)),
                                 "train.weighted"]) < 10^(-6)) &&
                loglik[ii - 1, "train.weighted"] > -Inf)) { 
          converged.tmp <- TRUE
        }
      }
    }
    
    if (!more.sparse.tmp){
      # if new sparser weights are actually not sparse, then merge clusters
      do.merge.tmp <- TRUE
    }
    if (control$trace > 0) {    
      cat("Solution more sparse:", more.sparse.tmp, "\n")
      cat("EM converged temporarily:", converged.tmp, "\n")
      cat("Should states be merged:", do.merge.tmp, "\n")
    }
    loglik[ii, "train.weighted"] <-
      compute_LICORS_loglik(pdfs$FLC$train, weight.matrices$train$tmp,
                            lambda = control$lambda)
    # TODO: make this state conversion in one step; not calling two functions
    loglik[ii, "train"] <-
        compute_LICORS_loglik(pdfs$FLC$train, 
                              states2weight_matrix(weight_matrix2states(weight.matrices$train$tmp),
                                                   num.states.total =
                                                     num.states.tmp),
                              lambda = control$lambda )
    
    penalty[ii,"L2"] <- compute_mixture_penalty(weight.matrices$train$tmp,
                                                  "Lq", q = 2)
    penalty[ii,"entropy"] <- compute_mixture_penalty(weight.matrices$train$tmp, 
                                                     "entropy",
                                                     base = "num.states")
    
    penalty[ii, "nec"] <- 
      log(num.states.tmp) * penalty[ii, "entropy"] /
        (loglik[ii, "train"] - loglik.one.state)
    penalty[ii, "nec.weighted"] <- log(num.states.tmp) *
        penalty[ii, "entropy"] /
            (loglik[ii, "train.weighted"] - loglik.one.state)
    
    if (do.merge.tmp) {
        AA <- estimate_state_adj_matrix(pdfs.FLC = pdfs$FLC$train, 
                                        alpha = NULL, 
                                        distance = 
                                          function(f, g) {
                                            return(mean(abs(f - g))) })
        
        diag(AA) <- 0 # for norm and testing (similarity)
        # if EM converged and no merging is possible, then stop iterations
        if (all(AA < control$alpha)){
          #if (TRUE){
          stop.iterations <- TRUE
        } else {
          merged <- TRUE
          merging.iter <- c(merging.iter, ii)
          weight.matrices$train$tmp <- merge_states(which(AA == max(AA),
                                                          arr.ind = TRUE)[1:2],
                                                    weight.matrices$train$tmp)
          
          converged.tmp <- FALSE
          if (control$trace > 0){
            cat("Merged two states into one after convergence with",
                num.states.tmp ," states. \n")
          }
        }
    }
    
    # keep only those states that have an effective sample size of at least 1
    effective.sample.sizes <- colSums(weight.matrices$train$tmp)
    too.small.states <- effective.sample.sizes < 1
    if (any(too.small.states)) {
      size.zero.states <- which(too.small.states)
      non.zero.size.state <- which(!too.small.states)
      weight.matrices$train$tmp <- 
        merge_states(c(non.zero.size.state[1], size.zero.states),
                     weight.matrices$train$tmp)
      merged <- TRUE
      converged.tmp <- FALSE
      
      if (control$trace > 0){
        cat("Merged states into one since they had 0 sample size. \n")
      }      
      if (is.null(dim(weight.matrices$train$tmp))) {
        stop.iterations <- TRUE
      }
      if (stop.iterations){
        num.iter["final"] <- ii
        overall.converged <- TRUE
        if (control$trace > 0) {
          cat("***************************************************** \n")
          cat("EM algorithm converged to (local) optimum at iteration",
              num.iter["final"],".\n")
          cat("***************************************************** \n \n")
        }
        break
      }
    }
    
    pdfs$PLC$train <- estimate_LC_pdfs(LCs = train.data$PLC, 
                                       weight.matrix =
                                         weight.matrices$train$tmp,
                                       method = control$estimation.method[["PLC"]])
    
    pdfs$FLC$train <- estimate_LC_pdfs(LCs = train.data$FLC, 
                                       weight.matrix =
                                         weight.matrices$train$tmp,
                                       method = control$estimation.method[["FLC"]])    

    loglik[ii, "train.weighted.update"] <-
      compute_LICORS_loglik(pdfs$FLC$train, weight.matrices$train$tmp, 
                            lambda = control$lambda)
    num.states.tmp <- ncol(weight.matrices$train$tmp)
    
    # evaluate in-sample MSE
    prediction$FLC$train <-
      predict_FLC_given_PLC(train = list(data = list(FLC = train.data$FLC),
                                         weight.matrix = weight.matrices$train$tmp),
                            test = list(weight.matrix = weight.matrices$train$tmp),
                            type = "mean")

    prediction$FLC$train.weighted <-
      predict_FLC_given_PLC(train = list(data = train.data,
                                         weight.matrix = weight.matrices$train$tmp),
                            test = list(weight.matrix = weight.matrices$train$tmp),
                            type = "weighted.mean")

    MSE[ii, "train.weighted"] <- control$loss(train.data$FLC, prediction$FLC$train.weighted)
    MSE[ii, "train"] <- control$loss(train.data$FLC, prediction$FLC$train)
    

    # evaluate out-of-sample MSE
    pdfs$PLC$test <- estimate_LC_pdfs(LCs = train.data$PLC, 
                                      weight.matrix = weight.matrices$train$tmp,
                                      method = control$estimation.method[["PLC"]],
                                      eval.LCs = test.data$PLC)
    weight.matrices$test$tmp <-
        estimate_state_probs(weight.matrix = weight.matrices$train$tmp,
                             pdfs = list(PLC = pdfs$PLC$test, FLC = NULL),
                             num.states = ncol(pdfs$PLC))
    state.vector$test <- weight_matrix2states(weight.matrices$test$tmp)

    prediction$FLC$test.weighted <-
      predict_FLC_given_PLC(train = list(data = train.data,
                                         weight.matrix = weight.matrices$train$tmp),
                            test = list(weight.matrix = weight.matrices$test$tmp),
                            type = "weighted.mean")
  
    prediction$FLC$test <-
      predict_FLC_given_PLC(train = list(data = train.data,
                                         weight.matrix = weight.matrices$train$tmp),
                            test = list(weight.matrix = weight.matrices$test$tmp),
                            type = "mean")

    MSE[ii, "test.weighted"] <- control$loss(test.data$FLC, prediction$FLC$test.weighted)
    MSE[ii, "test"] <- control$loss(test.data$FLC, prediction$FLC$test)
    pdfs$FLC$test <- estimate_LC_pdfs(LCs = train.data$FLC, 
                                      weight.matrix =
                                        weight.matrices$train$tmp,
                                      method = control$estimation.method[["FLC"]],
                                      eval.LCs = test.data$FLC)
    
    num.states["current"] <- ncol(pdfs$FLC$test)
    loglik[ii, "test.weighted"] <-
      compute_LICORS_loglik(pdfs$FLC$test, weight.matrices$test$tmp, 
                            lambda = control$lambda)
    loglik[ii, "test"] <-
      compute_LICORS_loglik(pdfs$FLC$test, 
                            states2weight_matrix(state.vector$test, 
                                                 num.states.total = num.states["current"]), 
                            lambda = control$lambda)
    if (MSE[ii, "test.weighted"] <= min(MSE[-ii, "test.weighted"], na.rm = TRUE)) {
      weight.matrices$train$opt <- weight.matrices$train$tmp
      weight.matrices$test$opt <- weight.matrices$test$tmp
      prediction$FLC$opt <- prediction$FLC$test
      prediction$FLC$opt.weighted <- prediction$FLC$test.weighted
      num.iter["opt"] <- ii
    }

    if (control$trace > 0) {
      cat("\n Finished step", ii, "with", num.states.tmp, "states. \n")
      cat("Loglik:", loglik[ii, "train.weighted"], "\n")
      cat("Out-of-sample loglik:", loglik[ii, "test.weighted"], "\n" )
      cat("Penalty:", round(penalty[ii, ]*100,1), "% \n")
      if (!merged){
        cat("Maximum difference in weights:", round(max.weight.diff, 4), "\n")
      } else {
  	    cat("Weights have been merged \n")
      }
      cat("In-sample MSE:", MSE[ii, "train"], "\n" )
      cat("Weighted In-sample MSE:", MSE[ii, "train.weighted"], "\n" )
      cat("Out-of-sample MSE:", MSE[ii, "test"], "\n" )
      cat("Weighted Out-of-sample MSE:", MSE[ii, "test.weighted"], "\n \n" )
      cat("In-sample data vs fit correlation:",
          round( 100 * cor(train.data$FLC, prediction$FLC$train.weighted), 1),
          "%\n")
      cat("Out-of-sample data vs fit correlation:",
          round( 100 * cor(test.data$FLC, prediction$FLC$test.weighted), 1),
          "%\n \n")
    }
    
    weight.matrices$train$final <- weight.matrices$train$tmp
   
    if (ncol(weight.matrices$train$tmp) < 3) {
      stop.iterations <- TRUE
    }

    if (stop.iterations){
      num.iter["final"] <- ii
      overall.converged <- TRUE
      if (control$trace > 0) {
        cat("***************************************************** \n")
        cat("EM algorithm converged to (local) optimum at iteration",
            num.iter["final"],".\n")
        cat("***************************************************** \n \n")
      }
      break
    }
    
  }
  
  if (!stop.iterations) {
    if (control$trace > 0) {
  	  cat("***************************************************** \n")
  	  cat("Finished all iterations without convergence. \n")
  	  cat("***************************************************** \n")
    }
  }
  
  weight.matrices$test$final <- weight.matrices$test$tmp
  prediction$states$test$final <-
    weight_matrix2states(weight.matrices$test$final)
  prediction$states$test$opt <-
    weight_matrix2states(weight.matrices$test$opt)

  num.states[c("opt", "final")] <- c(ncol(weight.matrices$test$opt),
                                     ncol(weight.matrices$test$final))
  num.iter["final"] <- ii
  
  out <- list(LCs = LCs,
              train.index = train.index,
              loglik.trace = loglik[seq_len(num.iter["final"]), ],
              loglik = loglik[c(1, num.iter["opt"], num.iter["final"]), ])
  rownames(out$loglik) <- c("init", "opt", "final")

  out <- c(out,
           list(MSE.trace = MSE[seq_len(num.iter["final"]), ],
                MSE = MSE[num.iter["opt"], ]))
  out$weight.matrices <- weight.matrices
  
  out$states <- rep(NA, NN)
  out$states[train.index] <- weight_matrix2states(weight.matrices$train$opt)
  out$states[-train.index] <- prediction$states$test$opt

  out$conditional.state.probs <- list(opt = NULL, final = NULL)
   out$conditional.state.probs$final <- 
    Matrix(0, ncol = num.states["final"], nrow = num.samples["all"],
           sparse = TRUE)
  # copy format to optimal matrix
  out$conditional.state.probs$opt <- 
    Matrix(0, ncol = num.states["opt"], nrow = num.samples["all"],
           sparse = TRUE)

  out$conditional.state.probs$final[train.index, ] <- weight.matrices$train$final
  out$conditional.state.probs$final[-train.index, ] <- weight.matrices$test$final
  out$conditional.state.probs$opt[-train.index, ] <- weight.matrices$test$opt
  out$conditional.state.probs$opt[train.index, ] <- weight.matrices$train$opt

  out$marginal.state.probs <- lapply(out$conditinal.state.probs,
                                     estimate_state_probs)
  out$penalty <- penalty
  out$loglik.one.state <- loglik.one.state

  out$num.states <- num.states["current" != names(num.states)]
  out$num.iter <- num.iter

  out$control <- control
  out$merging.iter <- merging.iter
  out$num.merges <- length(out$merging.iter)
  
  out$converged <- overall.converged
  out$dim <- LCs$dim
  class(out) <- c("LICORS", "mixed_LICORS")
  invisible(out)
}
