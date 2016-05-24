#' Variable Importance Measure for A FSM Model
#'
#' \code{var_imp} calculates the importance of the covariates of the model.
#'
#' Takes the state matrix and action vector from an already evolved model and
#' the fitness function and data used to evolve the model (or this could be test
#' data), flips the values of each of the elements in the state matrix and
#' measures the change in fitness (prediction of data) relative to the original
#' model. Then these changes are summed across columns to provide the importance
#' of each of the columns of the state matrix.
#'
#' @param state_mat Numeric matrix with rows as states and columns as
#'   predictors.
#' @param action_vec Numeric vector indicating what action to take for each
#'   state.
#' @param data Data frame that has "period" and "outcome" columns and rest of
#'   cols are predictors, ranging from one to three predictors. All of the (3-5
#'   columns) should be named.
#' @param outcome Numeric vector same length as the number of rows as data.
#' @param period Numeric vector same length as the number of rows as data.
#' @param measure Optional length one character vector that is either:
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
#'
#' @return Numeric vector the same length as the number of columns of the
#'   provided state matrix (the number of predictors in the model) with relative
#'   importance scores for each predictor.
#'
#' @export

var_imp <- function(state_mat, action_vec, data, outcome, period, measure){

  counter <- 1
  indices <- as.list(rep(NA, length(as.vector(state_mat))))
  for (i in seq(nrow(state_mat))){
    for (j in seq(ncol(state_mat))){
      indices[[counter]] <- c(i, j)
      counter <- counter + 1
    }
  }

  if (nrow(state_mat)*ncol(state_mat) != length(indices)){
    stop("Error in var_imp computation: 
         Numbers of elements of state matrix does not equal length of list to hold indices for each of those elements.")
  }

  fitness_mat <- state_mat

  results1 <- fitnessCPP(action_vec, state_mat, data, period)
  if (anyNA(results1) | length(results1)==0){
    stop("Error in var_imp computation: 
         Results from initial fitness evaluation have missing values or are wrong length.")
  }
  
  results1 <- performance(results = results1, outcome = outcome, measure = measure)

  for (i in seq(length(indices))) {
    state_mat_flipped <- state_mat
    if(nrow(state_mat)==2){
      # If 2 rows in state matrix then just flip each value:
      state_mat_flipped[indices[[i]][1],
                        indices[[i]][2]] <- ifelse(state_mat[indices[[i]][1],
                                                             indices[[i]][2]]==1, 2, 1)
    } else {
      # Sample from all possible values of states besides the one you're in now:
      state_mat_flipped[indices[[i]][1],
                        indices[[i]][2]] <- sample(seq(nrow(state_mat))[-i], 1)
    }
    
    results2 <- fitnessCPP(action_vec, state_mat_flipped, data, period)
    
    if (anyNA(results2) | length(results2)==0){
      stop("Error in var_imp computation: 
           Results from subsequent fitness evaluation have missing values.")
    }

    results2 <- performance(results = results2, outcome = outcome, measure = measure)

    dif <- results1 - results2

    fitness_mat[indices[[i]][1],
                indices[[i]][2]] <- dif
  }

  varImp <- colSums(fitness_mat) # same as: as.vector(apply(fitness_mat, MARGIN=2, sum))
  varImp <- (varImp/ifelse(max(varImp)==0, 0.001, max(varImp)))*100 # make the best be 100
  varImp <- ifelse(varImp < 0, 0, varImp)
  names(varImp) <- colnames(state_mat)
  varImp
}
