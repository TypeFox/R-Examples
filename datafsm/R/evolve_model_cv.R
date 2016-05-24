
################################################################################
add_interact_num <- function(d){
        game <- rep(NA, nrow(d))
        game[1] <- 1
        for (i in 2:nrow(d)){
                game[i] <- ifelse(d$period[i]==1, game[i-1] + 1, game[i-1])
        }
        game
}


#' Estimate Optimal Number of States of a Finite-state Machine Model
#' 
#' \code{evolve_model_cv} calls \code{evolve_model} with varied numbers of 
#' states and compares their performance with cross-validation.
#' \code{evolve_model} relies on the \strong{GA} package for genetic algorithm
#' optimization.
#' 
#' @param data data.frame that has "period" and "outcome" columns and rest of 
#'   columns are predictors, ranging from one to three predictors. All of the 
#'   (3-5 columns) should be named.
#' @param measure Optional character vector that is either: "accuracy", "sens", 
#'   "spec", or "ppv". This specifies what measure of predictive performance to 
#'   use for training and evaluating the model. The default measure is 
#'   "accuracy". However, accuracy can be a problematic measure when the classes
#'   are imbalanced in the samples, i.e. if a class the model is trying to 
#'   predict is very rare. Alternatives to accuracy are available that 
#'   illuminate different aspects of predictive power. Sensitivity answers the 
#'   question, `` given that a result is truly an event, what is the probability
#'   that the model will predict an event?'' Specificity answers the question, 
#'   ``given that a result is truly not an event, what is the probability that 
#'   the model will predict a negative?'' Positive predictive value answers, 
#'   ``what is the percent of predicted positives that are actually positive?''
#' @param k Numeric vector with the number of folds for k-fold cross-validation.
#' @param max_states Numeric vector with the maximum number of states to test.
#' @param actions Numeric vector with the number of actions. If not provided, 
#'   then actions will be set as the number of unique values in the outcome 
#'   vector.
#' @param seed Numeric vector length one.
#' @param popSize Numeric vector length one specifying the size of the GA 
#'   population. A larger number will increase the probability of finding a very
#'   good solution but will also increase the computation time. This is passed 
#'   to the GA::ga() function of the \strong{GA} package.
#' @param pcrossover Numeric vector length one specifying probability of 
#'   crossover for GA. This is passed to the GA::ga() function of the 
#'   \strong{GA} package.
#' @param pmutation Numeric vector length one specifying probability of mutation
#'   for GA. This is passed to the GA::ga() function of the \strong{GA} package.
#' @param maxiter Numeric vector length one specifying max number of iterations 
#'   for stopping the GA evolution. A larger number will increase the 
#'   probability of finding a very good solution but will also increase the 
#'   computation time. This is passed to the GA::ga() function of the 
#'   \strong{GA} package.
#' @param run Numeric vector length one specifying max number of consecutive 
#'   iterations without improvement in best fitness score for stopping the GA 
#'   evolution. A larger number will increase the probability of finding a very 
#'   good solution but will also increase the computation time. This is passed 
#'   to the GA::ga() function of the \strong{GA} package.
#' @param parallel Logical vector length one. For running the GA evolution in 
#'   parallel. Depending on the number of cores registered and the memory on 
#'   your machine, this can make the process much faster, but only works for 
#'   Unix-based machines that can fork the processes.
#' @param verbose Optional logical vector length one specifying whether helpful 
#'   messages should be displayed on the user's console or not.
#'   
#' @return Returns the number of states that maximizes the \code{measure}, e.g.
#'   accuracy.
#'   
#' @references Luca Scrucca (2013). GA: A Package for Genetic Algorithms in R. 
#'   Journal of Statistical Software, 53(4), 1-37. URL 
#'   http://www.jstatsoft.org/v53/i04/.
#'   
#'   Hastie, T., R. Tibshirani, and J. Friedman. (2009). The Elements of
#'   Statistical Learning: Data Mining, Inference, and Prediction, Second
#'   Edition. 2nd ed. New York, NY: Springer.
#'   
#' @export

################################################################################
evolve_model_cv <- function(data,
                            measure,
                            k,
                            actions,
                            max_states,
                            seed,
                            popSize, pcrossover, pmutation, maxiter, run,
                            parallel,
                            verbose) {
  
  interacts <- add_interact_num(data)
  
  mat <- matrix(NA, nrow = max_states, ncol = k)
  
  for(s in seq(from = 2, to = max_states, by = 1)){
    # divide interacts into k folds
    group_folds <- caret::createFolds(y = unique(interacts), k = k, list = FALSE)
    
    if(length(group_folds) != length(unique(interacts))) stop("Assignment of groups to folds didnt work: length(group_folds) != length(unique(interacts)).")
    if(length(unique(group_folds)) != k) stop("Assignment of groups to folds didnt work: length(unique(group_folds)) != k.")
    # if(verbose) cat("Group folds:", group_folds ,"\n")
    
    # create a vector same length as data with assignments of each row to a fold:
    fold_ass <- rep(NA, nrow(data))
    for (i in seq(nrow(data))) fold_ass[i] <- group_folds[interacts[i]]
    if(length(fold_ass) != nrow(data))
      stop("Creating a vector same length as data with assignments of each row to a fold didnt work: length(fold_ass) != nrow(data).")
    if(length(unique(fold_ass)) != length(unique(group_folds)))
      stop("Creating a vector same length as data with assignments of each row to a fold didnt work: length(unique(fold_ass)) != length(unique(group_folds)).")
    
    # In the fth fold, the elements of folds that equal f are in the test set, and the remainder are in the training set.
    for(f in seq(k)){
      training <- fold_ass == f
      if(class(training) != "logical") stop("Training index not logical vector.")
      if(verbose) cat("\nCross-validated testing with states set to", s, "\n\n")
      mat[s, f] <- evolve_model(data[training, ], data[!training, ],
                                drop_nzv = FALSE,
                                measure = measure,
                                states = s, cv = FALSE, seed = seed,
                                popSize = popSize, pcrossover =pcrossover,
                                pmutation = pmutation, maxiter = maxiter, run = run,
                                parallel = parallel)@predictive
      if(verbose) cat("\nCross-validated value of", measure, "is", mat[s, f], ".\n\n")
    }
  }
  results <- apply(mat[seq(from = 2, to = max_states, by = 1), ], 1, mean) # mean predictive accuracy for each number of states across all k folds (columns)
  min(which(results==max(results))+1) # na.omit dropped the first row of mat bc we started at states==2
  # the number of states that maximizes accuracy obtained from index with highest value, but add one because
  # first position in vector corresponds to states==2
  # BREAKS TIES BY CHOOSING THE SMALLER NUMBER OF STATES (min(...))
}

# data = cdata; k=2; actions=2; max_states=4; seed=1; popSize = 75; pcrossover = 0.8; pmutation = 0.1;
# maxiter = 55; run = 25; parallel = FALSE
