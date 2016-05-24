sofia <- function(x, ...) {
  UseMethod("sofia")
}

sofia.formula <- function(x, data
  , random_seed = floor(runif(1, 1, 65535))
  , lambda = 0.1 
  , iterations = 100000
  , learner_type = c("pegasos", "sgd-svm", "passive-aggressive", "margin-perceptron", "romma", "logreg-pegasos") 
  , eta_type = c("pegasos", "basic", "constant") 
  , loop_type = c("stochastic", "balanced-stochastic", "rank", "roc", "query-norm-rank", "combined-ranking", "combined-roc")
  , rank_step_probability = 0.5
  , passive_aggressive_c = 10000000.0
  , passive_aggressive_lambda = 0
  , perceptron_margin_size = 1.0
  , training_objective = FALSE
  , hash_mask_bits = 0
  , verbose = FALSE
  , reserve = 0 
  , ...
) {



  ### need to replace with as.svmlight to eliminate duplicate code
  ### I think we should set no_bias_term permanately to FALSE and just have 
  ### the user provide the data with or without a column of ones

  ### for consistency the first argument needs to be x...

  learner_type <- match.arg(learner_type)
  loop_type    <- match.arg(loop_type)
  eta_type     <- match.arg(eta_type)

  if(class(x) != "formula")
    stop("x must be a formula")
  if(!is.data.frame(data))
    stop("data must be a dataframe")
  if(!(is.numeric(random_seed) && length(random_seed) == 1))
    stop("random_seed must be a numeric scalar")
  if(!(is.numeric(lambda) && length(lambda) == 1))
    stop("lambda must be a numeric scalar")
  if(!(learner_type %in% c("pegasos", "sgd-svm", "passive-aggressive", "margin-perceptron", "romma", "logreg-pegasos")))
    stop(sprintf("learner_type: %s not supported", learner_type))
  if(!(eta_type %in% c("pegasos", "basic", "constant")))
    stop(sprintf("eta_type: %s not supported", eta_type))
  if(!(loop_type %in% c("stochastic", "balanced-stochastic", "rank", "roc", "query-norm-rank", "combined-ranking", "combined-roc")))
    stop(sprintf("loop_type: %s not supported", loop_type))
  if(!(rank_step_probability >= 0 && rank_step_probability <= 1))
    stop("rank step probability must be between 0 and 1")
  if(!(is.numeric(passive_aggressive_c) && length(passive_aggressive_c)==1))
    stop("passive_aggressive_c must be a numeric scalar")
  if(!(is.numeric(passive_aggressive_lambda) && length(passive_aggressive_lambda)==1))
    stop("passive_aggressive_lambda must be a numeric scalar")
  if(!(is.numeric(perceptron_margin_size) && length(perceptron_margin_size)==1))
    stop("perceptron_margin_size must be a numeric scalar")
  if(!(is.logical(training_objective)))
    stop("training_objective must be 'TRUE' or 'FALSE'")
  if(!(is.numeric(hash_mask_bits) && length(hash_mask_bits)==1))
    stop("hash_mask_bits must be a numeric scalar")

  ### does verbose do anything??

  parsed <- parse_formula(x, data)

  dimensionality <- ncol(parsed$data)+1
  
  sofia.model <- sofia.fit(parsed$data, parsed$labels, random_seed, lambda, iterations, learner_type, eta_type, loop_type, rank_step_probability 
    , passive_aggressive_c, passive_aggressive_lambda, perceptron_margin_size, training_objective
    , parsed$no_bias_term, dimensionality, hash_mask_bits, verbose, reserve
  )
  
  sofia.model$formula <- x 
                  
  return(sofia.model)
                  
}

sofia.character <- function(x
  , random_seed = floor(runif(1, 1, 65535))
  , lambda = 0.1 
  , iterations = 100000
  , learner_type = c("pegasos", "sgd-svm", "passive-aggressive", "margin-perceptron", "romma", "logreg-pegasos") 
  , eta_type = c("pegasos", "basic", "constant") 
  , loop_type = c("stochastic", "balanced-stochastic", "rank", "roc", "query-norm-rank", "combined-ranking", "combined-roc")
	, rank_step_probability = 0.5
  , passive_aggressive_c = 10000000.0
  , passive_aggressive_lambda = 0
  , perceptron_margin_size = 1.0
  , training_objective = FALSE
  , no_bias_term = FALSE
  , dimensionality = 150000 
  , hash_mask_bits = 0
  , verbose = FALSE
  , buffer_mb = 40, ...) 
{
  sofia_facade <- new(RSofiaFacade)
  
  learner_type <- match.arg(learner_type)
  loop_type    <- match.arg(loop_type)
  eta_type     <- match.arg(eta_type)
  
  sofia_resultset <- sofia_facade$train_filename(x
    , random_seed 
    , lambda 
    , iterations
    , learner_type 
    , eta_type
    , loop_type
  	, rank_step_probability
    , passive_aggressive_c
    , passive_aggressive_lambda
    , perceptron_margin_size
    , training_objective
    , dimensionality
  	, hash_mask_bits
  	, no_bias_term
    , verbose
    , buffer_mb
  )

  weights        <- sofia_resultset$weights
  names(weights) <- c("(Offset)", seq_len(length(weights))[-length(weights)])
                                        
  training_time <- sofia_resultset$training_time 
  io_time       <- sofia_resultset$io_time
 
  obj <- list(
    par = list(random_seed=random_seed
                , lambda=lambda
                , iterations=iterations
                , learner_type=learner_type
                , eta_type=eta_type
                , loop_type=loop_type
                , rank_step_probability=rank_step_probability
                , passive_aggressive_c=passive_aggressive_c
                , passive_aggressive_lambda=passive_aggressive_lambda
                , perceptron_margin_size=perceptron_margin_size
                , training_objective=training_objective
                , dimensionality=dimensionality
              	, hash_mask_bits=hash_mask_bits
              	, no_bias_term=no_bias_term
                ),
    weights = weights
    , training_time = training_time
    , io_time       = io_time
  )
 
  class(obj) <- "sofia"
  
  return (obj)

}

sofia.fit <- function(x, y
  , random_seed = floor(runif(1, 1, 65535))
  , lambda = 0.1 
  , iterations = 100000
  , learner_type = "pegasos"
  , eta_type = "pegasos"
	, loop_type = "stochastic"
	, rank_step_probability = 0.5
  , passive_aggressive_c = 10000000.0
  , passive_aggressive_lambda = 0
  , perceptron_margin_size = 1.0
  , training_objective = FALSE
  , no_bias_term = FALSE
  , dimensionality = ncol(x) + 1
  , hash_mask_bits = 0
  , verbose = FALSE
  , reserve = 0 
  , ...
) {
  
  ###
  # break on bad parameter
  ###

  sofia_facade <- new(RSofiaFacade)
  
  sofia_resultset <- sofia_facade$train_fit(x, y
    , random_seed 
    , lambda 
    , iterations
    , learner_type 
    , eta_type
    , loop_type
  	, rank_step_probability
    , passive_aggressive_c
    , passive_aggressive_lambda
    , perceptron_margin_size
    , training_objective
    , dimensionality
  	, hash_mask_bits
  	, no_bias_term
    , verbose
    , reserve
  )

  weights        <- sofia_resultset$weights

  if(is.null(colnames(x))) {
    colnames_ <- as.character(seq_len(ncol(x)))
  } else {
    colnames_ <- colnames(x)
  }                 
 
  names(weights) <- c("(Offset)", colnames_) 
                                        
  training_time <- sofia_resultset$training_time 
  io_time       <- sofia_resultset$io_time
 
  obj <- list(
    par = list(random_seed=random_seed
                , lambda=lambda
                , iterations=iterations
                , learner_type=learner_type
                , eta_type=eta_type
                , loop_type=loop_type
                , rank_step_probability=rank_step_probability
                , passive_aggressive_c=passive_aggressive_c
                , passive_aggressive_lambda=passive_aggressive_lambda
                , perceptron_margin_size=perceptron_margin_size
                , training_objective=training_objective
                , dimensionality=dimensionality
              	, hash_mask_bits=hash_mask_bits
              	, no_bias_term=no_bias_term
                ),
    weights = weights
    , training_time = training_time
    , io_time = io_time
  )
 
  class(obj) <- "sofia"
  
  return (obj)

}
