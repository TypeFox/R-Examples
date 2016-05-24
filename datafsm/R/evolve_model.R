performance <- function(results, outcome, measure){
  switch(measure,
         accuracy = sum(ifelse( results == outcome , 1 , 0)) / length(results),
         sens = caret::sensitivity(data = factor(results), reference = factor(outcome)),
         spec = caret::specificity(data = factor(results), reference = factor(outcome)),
         ppv = caret::posPredValue(data = factor(results), reference = factor(outcome)))
}

#'Use a Genetic Algorithm to Estimate a Finite-state Machine Model
#'
#'\code{evolve_model} uses a genetic algorithm to estimate a finite-state 
#'machine model, primarily for understanding and predicting decision-making.
#'
#'This is the main function of the \strong{datafsm} package. It relies on the 
#'\strong{GA} package for genetic algorithm optimization. \code{evolve_model} 
#'takes data on predictors and data on the outcome. It automatically creates a 
#'fitness function that takes the data, an action vector \code{evolve_model} 
#'generates, and a state matrix \code{evolve_model} generates as input and 
#'returns numeric vector of the same length as the \code{outcome}. 
#'\code{evolve_model} then computes a fitness score for that potential solution 
#'FSM by comparing it to the provided \code{outcome}. This is repeated for every
#'FSM in the population and then the probability of selection for the next 
#'generation is proportional to the fitness scores. The default is also for the 
#'function to call itself recursively while varying the number of states inside 
#'a cross-validation loop in order to estimate the optimal number of states.
#'
#'If parallel is set to TRUE, then these evaluations are distributed across the 
#'available processors of the computer using the \strong{doParallel} package, 
#'otherwise, the evalulations of fitness are conducted sequentially. Because 
#'this fitness function that \code{evolve_model} creates must loop through all 
#'the data everytime it is evaluated and we need to evaluate many possible 
#'solution FSMs, the fitness function is implemented in C++ so it is very fast.
#'
#'\code{evolve_model} uses a stochastic meta-heuristic optimization routine to 
#'estimate the parameters that define a FSM model. Generalized simulated 
#'annealing, or tabu search could work, but they are more difficult to 
#'parallelize. The current version uses the \strong{GA} package's genetic 
#'algorithm because GAs perform well in rugged search spaces to solve integer 
#'optimization problems, are a natural complement to our binary string 
#'representation of FSMs, and are easily parallelized.
#'
#'This function evolves the models on training data and then, if a test set is 
#'provided, uses the best solution to make predictions on test data. Finally, 
#'the function returns the GA object and the decoded version of the best string 
#'in the population. See \linkS4class{ga_fsm} for the details of the slots 
#'(objects) that this type of object will have.
#'
#'@param data data.frame that has columns named "period" and "outcome" (period 
#'  is the time period that the outcome action was taken), and the rest of the 
#'  columns are predictors, ranging from one to three predictors. All of the 
#'  (3-5 columns) should be named. The period and outcome columns should be 
#'  integer vectors and the columns with the predictor variable data should be 
#'  logical vectors (\code{TRUE, FALSE}). If the predictor variable data is not 
#'  logical, it will coerced to logical with \code{base::as.logical()}.
#'@param test_data Optional data.frame that has "period" and "outcome" columns 
#'  and rest of columns are predictors, ranging from one to three predictors. 
#'  All of the (3-5 columns) should be named. Outcome variable is the decision 
#'  the decision-maker took for that period. This data.frame should be in the
#'  same format and have the same order of columns as the data.frame passed to
#'  the required \code{data} argument.
#'@param drop_nzv Optional logical vector length one specifying whether 
#'  predictors variables with variance in provided data near zero should be 
#'  dropped before model building. Default is \code{FALSE}. See
#'  \code{caret::nearZeroVar()}, which calls: \code{caret::nzv()}.
#'@param measure Optional length one character vector that is either:
#'  "accuracy", "sens", "spec", or "ppv". This specifies what measure of
#'  predictive performance to use for training and evaluating the model. The
#'  default measure is \code{"accuracy"}. However, accuracy can be a problematic
#'  measure when the classes are imbalanced in the samples, i.e. if a class the
#'  model is trying to predict is very rare. Alternatives to accuracy are
#'  available that illuminate different aspects of predictive power. Sensitivity
#'  answers the question, `` given that a result is truly an event, what is the
#'  probability that the model will predict an event?'' Specificity answers the
#'  question, ``given that a result is truly not an event, what is the
#'  probability that the model will predict a negative?'' Positive predictive
#'  value answers, ``what is the percent of predicted positives that are
#'  actually positive?''
#'@param states Optional numeric vector with the number of states. 
#'  If not provided, will be set to \code{max(data$outcome)}.
#'@param cv Optional logical vector length one for whether cross-validation 
#'  should be conducted on training data to select optimal number of states. 
#'  This can drastically increase computation time because if \code{TRUE}, it 
#'  will run evolve_model k*max_states times to estimate optimal value for 
#'  states. Ties are broken by choosing the smaller number of states. Default is
#'  \code{FALSE}.
#'@param max_states Optional numeric vector length one only relevant if 
#'  \code{cv==TRUE}. It specifies how up to how many states that 
#'  cross-validation should search through. 
#'  If not provided, will be set to \code{states + 1}.
#'@param k Optional numeric vector length one only relevant if cv==TRUE, 
#'  specifying number of folds for cross-validation.
#'@param actions Optional numeric vector with the number of actions. If not 
#'  provided, then actions will be set as the number of unique values in the 
#'  outcome vector.
#'@param seed Optional numeric vector length one.
#'@param popSize Optional numeric vector length one specifying the size of the 
#'  GA population. A larger number will increase the probability of finding a 
#'  very good solution but will also increase the computation time. This is 
#'  passed to the GA::ga() function of the \strong{GA} package.
#'@param pcrossover Optional numeric vector length one specifying probability of
#'  crossover for GA. This is passed to the GA::ga() function of the \strong{GA}
#'  package.
#'@param pmutation Optional numeric vector length one specifying probability of 
#'  mutation for GA. This is passed to the GA::ga() function of the \strong{GA} 
#'  package.
#'@param maxiter Optional numeric vector length one specifying max number of 
#'  iterations for stopping the GA evolution. A larger number will increase the 
#'  probability of finding a very good solution but will also increase the 
#'  computation time. This is passed to the GA::ga() function of the \strong{GA}
#'  package. \code{maxiter} is scaled by how many parameters are in the model: 
#'  \code{maxiter <- maxiter + ((maxiter*(nBits^2)) / maxiter)}.
#'@param run Optional numeric vector length one specifying max number of 
#'  consecutive iterations without improvement in best fitness score for 
#'  stopping the GA evolution. A larger number will increase the probability of 
#'  finding a very good solution but will also increase the computation time. 
#'  This is passed to the GA::ga() function of the \strong{GA} package.
#'@param parallel Optional logical vector length one. For running the GA 
#'  evolution in parallel. Depending on the number of cores registered and the 
#'  memory on your machine, this can make the process much faster, but only 
#'  works for Unix-based machines that can fork the processes.
#'@param priors Optional numeric matrix of solutions strings to be included in 
#'  the initialization. User needs to use a decoder function to translate prior 
#'  decision models into bits and then provide them. If this is not specified, 
#'  then random priors are automatically created.
#'@param verbose Optional logical vector length one specifying whether helpful 
#'  messages should be displayed on the user's console or not.
#'  
#'@return Returns an S4 object of class ga_fsm. See \linkS4class{ga_fsm} for the
#'  details of the slots (objects) that this type of object will have and for 
#'  information on the methods that can be used to summarize the calling and 
#'  execution of \code{evolve_model()}, including \code{summary}, 
#'  \code{print}, and \code{plot}. Timing measurement is in seconds.
#'  
#'@references Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R. 
#'  Journal of Statistical Software, 53(4), 1-37. URL 
#'  http://www.jstatsoft.org/v53/i04/.
#'  
#' @examples
#' # Create data:
#'cdata <- data.frame(period = rep(1:10, 1000),
#'                    outcome = rep(1:2, 5000),
#'                    my.decision1 = sample(1:0, 10000, TRUE),
#'                    other.decision1 = sample(1:0, 10000, TRUE))
#' (res <- evolve_model(cdata, cv=FALSE))
#' summary(res)
#' plot(res, action_label = c("C", "D"))
#' library(GA)
#' plot(estimation_details(res))
#'
#' # In scripts, it can makes sense to set parallel to
#' # 'as.logical(Sys.info()['sysname'] != 'Windows')'.
#'
#'@export

################################################################################
evolve_model <- function(data, test_data = NULL, drop_nzv = FALSE,
                         measure = c("accuracy", "sens", "spec", "ppv"),
                         states = NULL, cv = FALSE, max_states = NULL, k = 2,
                         actions = NULL,
                         seed = NULL,
                         popSize = 75, pcrossover = 0.8, pmutation = 0.1, maxiter = 50, run = 25,
                         parallel = FALSE,
                         priors = NULL,
                         verbose = TRUE) {
  
  start_time <- as.numeric(proc.time()[[3]])
  
  call <- match.call()
  
  msg <- "\n\n"
  
  measure <- match.arg(measure) # inside fitnessR(), use this arg to create measure of fitness
  
  ## GA-related errors:
  if (popSize < 10) warning("The population size is less than 10. Consider using a larger size.")
  if (maxiter < 1) stop("Error: The maximum number of iterations must be at least 1.")
  if (pcrossover < 0 || pcrossover > 1) stop("Error: Probability of crossover must be between 0 and 1.")
  if (pmutation < 0 || pmutation > 1) stop("Error: Probability of mutation must be between 0 and 1.")
  
  ## Parallel-related errors:
  if(!requireNamespace("doParallel", quietly = TRUE) & parallel == TRUE)
    stop(paste("You asked to run this in parallel, but you don't have package doParallel installed.",
               "run 'install.packages(''doParallel''); library(doParallel)', and then try this again."))
  
  ## Data-related errors:
  if (missing(data)) 
    stop(paste("You must supply data. At the very least, you can supply data.",
               "This should be a data.frame that has columns named 'period' and 'outcome' (period",
               "is the time period that the outcome action was taken), and the rest of the",
               "columns are predictors, ranging from one to three predictors. All of the",
               "(3-5 columns) should be named. The period and outcome columns should be",
               "integer vectors and the columns with the predictor variable data should be",
               "logical vectors."))
  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
    warning(paste("You did not supply a data.frame for 'data' argument of this function",
                  "so we converted it to one. To ensure this works right, run this again with a",
                  "data.frame that has columns named 'period' and 'outcome' (period",
                  "is the time period that the outcome action was taken), and the rest of the",
                  "columns are predictors, ranging from one to three predictors. All of the",
                  "(3-5 columns) should be named. The period and outcome columns should be",
                  "integer vectors and the columns with the predictor variable data should be",
                  "logical vectors."))
  }
  period <- data$period
  outcome <- data$outcome
  
  nzvs <- caret::nearZeroVar(data[ , -which(names(data) %in% c("period", "outcome")), drop=FALSE],
                             freqCut = 95/5, uniqueCut=10)
  if (length(nzvs) > 0){
    to_drop <- colnames(data)[-which(names(data) %in% c("period", "outcome"))[nzvs]]
    if(verbose) cat("We should be dropping", length(to_drop), "feature(s), which is (are):", to_drop, "\n")
    msg <- paste0(msg, "We should be dropping ", length(to_drop), " feature(s), which is (are): ", to_drop, "\n")
    
    if(drop_nzv){
      # just names in features[[k]] so we dont drop group, folds and training vars
      if(verbose) cat("Dropping", length(to_drop), "feature(s), which is (are):", to_drop)
      msg <- paste0(msg, "Dropping ", length(to_drop), " feature(s), which is (are): ", to_drop)
      
      data <- data[ , -which(names(data) %in% to_drop), drop=FALSE]
    }
  }
  
  if (nrow(data)!=length(outcome)) stop(paste("Error: The predictor variables and the",
                                              "outcome variable are not the same length."))
  if (anyNA(outcome)) stop("Error: There are missing values in the data.")
  if (length(outcome) == 0) stop("Error: The outcome is zero length.")
  if (missing(seed)) {
    seed <- floor(stats::runif(1, 1,101))
    if(verbose) cat(paste("We set a seed for you to make this reproducible. It is ", seed, ".",
                          " If you want the same results, next time you run this with the same settings,",
                          " also set the seed argument of this function to ", seed, ".\n", sep=""))
    msg <- paste0(msg, paste("We set a seed for you to make this reproducible. It is ", seed, ".",
                             " If you want the same results, next time you run this with the same settings,",
                             " also set the seed argument of this function to ", seed, ".\n", sep=""))
  }
  if (missing(actions)) {
    if(length(unique(outcome))==1){
      stop(paste("Error: There is only one unique value in the",
                 "outcome vector you supplied."))
    } else {
      actions <- length(unique(outcome))
    }
  } else {
    if (length(unique(outcome)) != actions) {warning(paste("The number of unique values in the",
                                                           "outcome vector you supplied does not",
                                                           "equal the value of actions you supplied.",
                                                           "The outcome vector should be a vector of",
                                                           "observed actions. We are going to use the",
                                                           "number of unique values in the outcome",
                                                           "vector you supplied as the value of actions."))
      actions <- length(unique(outcome))
    }
  }
  
  # So we are assured that the action vec will just need to be comprised of the possible
  # number of actions in the data:
  if (!identical(sort(as.integer(unique(outcome))), 
                 sort(as.integer(unique(seq(length(unique(outcome)))))),
                 ignore.environment = TRUE)){
    stop(paste("Error: The actions in the outcome column of the data are not the right values.",
               "There should be actions sequenced from 1 to however many actions that are feasible.",
               "E.g., if there are two feasible actions, then the outcome column should be comprised",
               "of only 1s and 2s, with at least one 1 and and at least one 2. If there are three feasible",
               "actions, the outcome column should be comprised of only 1s, 2s, and 3s, with at least one",
               "1 and, at least one 2, and at least one 3."))
  }
  
  inputs <- 2^(ncol(data[ , -which(names(data) %in% c("period", "outcome")), drop = FALSE]))
  
  if (is.null(states)) states <- max(data$outcome)
  if (is.null(max_states)) max_states <- states + 1

  if (cv) {
    try({
      states <- evolve_model_cv(data,
                                measure,
                                k,
                                actions,
                                max_states,
                                seed,
                                popSize, pcrossover, pmutation, maxiter, run,
                                parallel,
                                verbose)
      if(verbose) cat("Cross-validation found optimal number of states on training data to be ", states, ".\n\n", sep="")
      msg <- paste0(msg, "Cross-validation found optimal number of states on training data to be ", states, ".\n\n")
    })
    # wrapped this in try, so if it fails, we'll just use the default value of states, which is an arg to evolve_model()
  }
  
  # change any non-logical predictor variable vectors to logical
  data[ , -which(names(data) %in% c("period", "outcome"))] <-
    data.frame(lapply(data[ , -which(names(data) %in% c("period", "outcome"))],
                      function(x) {
                        if (class(x)!="logical") {
                          as.logical(x)
                        } else {
                          x
                        }}))
  # replace all NA's with 0 or 1 so these rows are not dropped
  # this works fine if the NAs are only for the first period play bc
  # then the predictor columns dont make a difference bc the FSM will initialize
  # with the same action regardless of the predictors at that time
  # but this would bias the results if NA's are occuring in predictors in other periods
  # so return an error for that:
  if (any(!stats::complete.cases(data) & !data$period == 1))
    stop(paste("Error: You have missing values in your training data somewhere other than the first period interactions.",
               "You can only have missing values for predictor columns, AND these must be in rows where period==1."))
  data[is.na(data)] <- TRUE
  
  names <- colnames(data[ , -which(names(data) %in% c("period", "outcome")), drop = FALSE])
  
  if (length(names)==1){
    form <- paste("outcome ~ 0 +", names, sep=" ")
    data <- stats::model.matrix(eval(parse(text=form)), data)
  } else {
    predictors <- paste(names, collapse=":")
    form <- paste("outcome ~ 0 +", predictors, sep=" ")
    data <- stats::model.matrix(eval(parse(text=form)), data)
  }
  
  if (length(names) > 3) stop(paste("Error: You have more than 3 predictors.",
                                    "Your model will be too complicated.",
                                    "Do some type of feature selection to choose less",
                                    "than 4 predictors and then use the data.frame",
                                    "with just those features next time."))
  
  if (ncol(data) != inputs)
    stop("Error: At least one of your predictor variables does not have exactly 2 levels.")
  
  cols <- colnames(data)
  # numeric vector same length as number of columns of the
  # state matrix (\code{state_mat}) with the action that each column of the
  # state matrix corresponds to the decision model taking in the previous
  # period. This is only relevant when the predictor variables of the FSM are
  # lagged outcomes that include the previous actions taken by that decision model.
  
  fitnessR <- function(s){ # Functions defined elsewhere in datafsm pkg: decode_action_vec, decode_state_mat, fitnessCPP
    action_vec <- decode_action_vec(s, states, inputs, actions)
    state_mat <- decode_state_mat(s, states, inputs, actions)
    results <- fitnessCPP(action_vec, state_mat, data, period)
    
    if (anyNA(results) | length(results)==0){
      stop("Error: Results from fitness evaluation have missing values.")
    }
    
    performance(results = results, outcome = outcome, 
                measure = measure)
  }
  
  warning_threshold <- 100
  
  valid_bs <- function(bs) {
    a <- decode_action_vec(bs, states, inputs, actions )
    sm <- decode_state_mat(bs, states, inputs, actions )
    all(a <= actions) && all(sm <= states)
  }
  
  valid_bsl <- function(x) {
    vbs <- function(i) valid_bs(x[i,])
    as.logical(lapply(1:nrow(x), vbs))
  }
  
  spCrossover <- function(object, parents, ...) {
    iter <- 0
    output <- NULL
    while(is.null(output)) {
      iter <- iter + 1
      output <- GA::gabin_spCrossover(object, parents, ...)
      children <- output$children
      if (! all(valid_bsl(children))) {
        if (iter > warning_threshold) {
          warning("Invalid crossover #", iter, '\n')
          print(children)
        }
        output <- NULL
      }
    }
    output
  }
  
  raMutation <- function(object, parent, ...) {
    iter <- 0
    output <- NULL
    while(is.null(output)) {
      iter <- iter + 1
      output <- GA::gabin_raMutation(object, parent, ...)
      if (! valid_bs(output)) {
        if (iter > warning_threshold) {
          cat("Invalid mutation #", iter, '\n')
          print(output)
        }
        output <- NULL
      }
    }
    output
  }
  
  poss.state.values <- seq(states) - 1
  b1 <- GA::decimal2binary(max(poss.state.values))
  l1 <- length(b1) #how many binary elements to represent one element of state matrix
  
  poss.action.values <- seq(actions) - 1
  b2 <- GA::decimal2binary(max(poss.action.values))
  l2 <- length(b2) #how many binary elements to represent one element of action matrix
  
  nBits <- (states*inputs*l1 + states*l2)
  
  build_priors <- function(popSize, nBits, states, inputs, actions) {
    priors <- matrix(nrow = popSize, ncol = nBits)
    for(i in 1:popSize) {
      av <- sample(actions,states,TRUE)
      sm <- matrix(sample(states,states * inputs,TRUE), nrow=states)
      priors[i,] <- build_bitstring(av,sm, actions)
    }
    #    prior_fitness <- unlist(lapply(1:popSize, function(i) fitnessR(priors[i,])))
    priors
  }
  
  if (missing(priors)) {
    priors <- build_priors(popSize, nBits, states, inputs, actions)
  } else {
    if (is.vector(priors)) {
      if (nBits > 1){
        priors <- matrix(priors, nrow = 1)
      } else  {
        priors <- matrix(priors, ncol = 1)
      }
    } else {
      priors <- as.matrix(priors)
    }
    if (nBits != ncol(priors)) stop(paste("Error: Priors do not match number of variables.",
                                          "Remember that you need to provide a decoded bitstring for the priors."))
  }
  
  # Scale convergence criteria by how many parameters are in the model:
  maxiter <- maxiter + ((maxiter*(nBits^2)) / maxiter)
  #run <- run + ((run*(nBits)) / run) # TODO: think about how to scale run
  
  if(verbose){
    Monitor <- function (object, digits = getOption("digits")) {
      Fit <- stats::na.exclude(object@fitness)
      cat(paste("Iter =", object@iter, " | Mean =", format(mean(Fit),
                                                           digits = digits), " | Best =", format(max(Fit),
                                                                                                 digits = digits),
                "\n"))
    }
  } else {
    Monitor <- function (object, digits = getOption("digits")) { }
  }
  
  GA <- GA::ga(type = "binary",
               fitness = fitnessR,
               nBits = nBits,
               crossover = spCrossover,
               mutation = raMutation,
               popSize = popSize,
               pcrossover = pcrossover,
               pmutation = pmutation,
               maxiter = maxiter,
               run = run,
               maxfitness = 1,
               parallel = parallel,
               suggestions = priors,
               monitor = Monitor,
               seed = seed)
  
  state_mat <- decode_state_mat(GA@solution[1, ],  states, inputs, actions)
  colnames(state_mat) <- cols
  
  action_vec <- decode_action_vec(GA@solution[1, ],  states, inputs, actions)
  
  if (missing(test_data)){
    predictive <- "No test data provided. Provide some to get more accurate estimation of generalization power."
  } else {
    
    test_period <- test_data$period
    test_outcome <- test_data$outcome
    
    if (!identical(sort(as.integer(unique(test_outcome))), 
                   sort(as.integer(unique(seq(length(unique(test_outcome)))))),
                   ignore.environment = TRUE)){
      warning(paste("Error: The actions in the outcome column of the test data are not the right values.",
                 "There should be actions sequenced from 1 to however many actions that are feasible.",
                 "E.g., if there are two feasible actions, then the outcome column should be comprised",
                 "of only 1s and 2s, with at least one 1 and and at least one 2. If there are three feasible",
                 "actions, the outcome column should be comprised of only 1s, 2s, and 3s, with at least one",
                 "1 and, at least one 2, and at least one 3."))
    }
    
    test_inputs <- 2^(ncol(test_data[ , -which(names(test_data) %in% c("period", "outcome")), drop = FALSE]))
    
    # change any non-logical predictor variable vectors to logical
    test_data[ , -which(names(test_data) %in% c("period", "outcome"))] <-
      data.frame(lapply(test_data[ , -which(names(test_data) %in% c("period", "outcome"))],
                        function(x) {
                          if (class(x)!="logical") {
                            as.logical(x)
                          } else {
                            x
                          }}))
    
    # replace all NA's with 0 or 1 so these rows are not dropped
    # this works fine if the NAs are only for the first period play bc
    # then the predictor columns dont make a difference bc the FSM will initialize
    # with the same action regardless of the predictors at that time
    # but this would bias the results if NA's are occuring in predictors in other periods
    # so return an error for that:
    if (any(!stats::complete.cases(test_data) & !test_data$period == 1))
      warning("Error: You have missing values in your test data somewhere other than the first period interactions. You can only have missing values for predictor columns, AND these must be in rows where period==1.")
    test_data[is.na(test_data)] <- TRUE
    
    names <- colnames(test_data[ , -which(names(test_data) %in% c("period", "outcome"))])
    
    if (length(names)==1){
      form <- paste("outcome ~ 0 +", names, sep=" ")
      test_data <- stats::model.matrix(eval(parse(text=form)), test_data)
    } else {
      predictors <- paste(names, collapse=":")
      form <- paste("outcome ~ 0 +", predictors, sep=" ")
      test_data <- stats::model.matrix(eval(parse(text=form)), test_data)
    }
    
    if (length(names) > 3) warning(paste("Error: You have more than 3 predictors in your test data.",
                                      "Your model will be too complicated.",
                                      "Do some type of feature selection to choose less",
                                      "than 4 predictors and then use the data.frame",
                                      "with just those features next time."))
    
    if (ncol(test_data) != test_inputs)
      warning(paste("Error: At least one of your predictor variables in your test data",
                 "does not have exactly 2 levels."))
    
    if(!all(colnames(test_data) == colnames(data)))
      warning("colnames(test_data) != colnames(data) and therefore test results will not be useful.")
    
    results <- fitnessCPP(action_vec, state_mat, test_data, test_period)
    if (anyNA(results) | length(results)==0){
      warning("Error: Results from fitness evaluation have missing values.")
    }
    
    predictive <- performance(results = results, outcome = test_outcome, 
                              measure = measure)
  }
  
  varImp <- var_imp(state_mat, action_vec, data, outcome, period, measure)
  
  methods::new("ga_fsm",
               call = call,
               actions = actions,
               states = states,
               GA = GA,
               state_mat = state_mat,
               action_vec = action_vec,
               predictive = predictive,
               varImp = varImp,
               timing = as.numeric(proc.time()[[3]]) - start_time,
               diagnostics = msg)
}
