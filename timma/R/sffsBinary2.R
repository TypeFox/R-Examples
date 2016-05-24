#' Model selection with filtered binary drug-target interaction data
#' 
#' A function to run sffs for model selection with filtered binary drug-target interaction data
#' 
#' @param profile_data drug-target interaction data which is a matrix with drugs as row indexes and targets
#' as column indexes.
#' @param sens a drug sensitivity vector.
#' @param sp an integer to specify the starting point for sequential forward floating search (sffs) search 
#' algorithm to navigate the target set space. By default, sp = 1.
#' @param max_k an integer to sepcify the maximum number of targets that can be selected by the sffs
#' algorithm. By default, max_k = 5. In practice it should not be over than 10 as the number of target combinations will increase exponentially.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#' @param new_initial_list a vector of the filtered targets indexes.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' 
#' @return A list containing the following components:
#' \item{timma}{a list contains: the predicted efficacy matrix, prediction error and predicted drug sensitivity}
#' \item{k_sel}{the indexes for selected targets}
#' @details The major difference between original and modified averaging method is the averaging methods for the case where the minimization and maximization rules are not simultaneously satisfied. 
#' For example, for a queried target set there are supersets but not subsets in the training data, the original algorithm will take the prediction from these supersets data using the minimization rule. 
#' However, the modified algorithm will further adjust the prediction using the average between such a prediction and 0. 
#' @author Liye He \email{liye.he@@helsinki.fi}
#' @examples
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' profile<-tyner_interaction_binary[,c(-1, -2, -5)]
#' num<-length(tyner_sensitivity[,1])
#' k_set<-rep(0, dim(profile)[2])
#' k_set[1]<-1
#' result<-sffsBinary2(profile, tyner_sensitivity[,1], new_initial_list = k_set, max_k=2)
#' }

sffsBinary2 <- function(profile_data, sens, sp = 1, max_k = 5, loo = TRUE, new_initial_list, verbosity = FALSE) {
  # sffs for binary timma model
  
  # get the number of the drugs
  drug_num <- nrow(profile_data)
  # get the number of the kinases
  kinase_num <- ncol(profile_data)
  
  # the criterion function J(X_k)
  J <- rep(0, kinase_num)
  # the set X_k
  X <- matrix(0, kinase_num, kinase_num)
  
  run <- 0
  step <- 0
  change <- 1
  global_best_err <- 100
  
  while (change) {
    run <- run + 1
    k_set <- new_initial_list
    
    while (sum(k_set) <= max_k) {
      err <- rep(Inf, kinase_num)
      step <- sum(k_set)
      profile_temp <- profile_data[, which(k_set == 1)]
      timma_temp <- timmaBinary(profile_temp, sens, loo)
      # the score for X_k, the lower score is better
      J[step] <- sum(timma_temp$error)
      # if k_set has been visited
      if (all(X[, step] == k_set)) {
        if(verbosity){
          cat("k_set has been visted!")
        }
        
        tmp <- which(J == min(J[which(J > 0)]))
        k_set <- X[, tmp[1]]
        change <- 0
        break
      } else {
        X[, step] <- k_set
      }
      
      if(verbosity){
        cat("Number of selected targets: ", step, "MAE = ", J[step], "\n")
      }
      
      
      if (step > 1 && J[step] == J[step - 1]) {
        # if no improvement
        new_initial_list <- rep(0, kinase_num)
        new_initial_list[ind1] <- 1
        if (J[step] < global_best_err) {
          global_best_err = J[step]
          change <- 1
          if(verbosity){
            cat("Best MAE:", global_best_err, "\n")
          }
          
        } else {
          change <- 0
        }
        break
      }
      
      # step 1: Inclusion
      space <- searchSpace(drug_num, k_set, profile_data, sens)
      for (i in 1:kinase_num) {
        if (k_set[i] == 0) {
          # including new kinase i
          err[i] <- 0
          err1 <- timmaSearchBinary(profile_data[, i], space, sens, loo)
          err[i] <- sum(err1)
        }
      }
      
      # select X_{k+1} the most significant feature with respect to X_k
      dummy <- min(err, na.rm = TRUE)
      ind1 <- which.min(err)
      k_set[ind1] <- 1
      if(verbosity){
        cat("Including target #", ind1, "MAE = ", dummy, "\n")
      }
      
      
      # step2: Conditional exclusion find the least significant feature in the set X_{k+1}
      err <- rep(Inf, kinase_num)
      for (i in 1:kinase_num) {
        k_temp <- k_set  # the current kinase set
        if (k_set[i] == 1) {
          # if removing the i-th kinase
          err[i] <- 0
          k_temp[i] <- 0
          profile_temp <- profile_data[, which(k_temp == 1)]  # which(k_temp==1)
          timma_temp <- timmaBinary(profile_temp, sens, loo)
          err[i] <- sum(timma_temp$error)
        }
      }
      worst_err <- min(err)
      ind_worst <- which.min(err)
      if (err[ind1] != worst_err && worst_err < J[step]) {
        # the new added kinase is not the least significant feature
        k_set[ind_worst] <- 0  # exclude the least significant feature
        if(verbosity){
          cat("Excluding target #", ind_worst, "MAE = ", worst_err, "\n")
        }
        
      }
      # step 3: Continuation of conditional exclusion
      while (sum(k_set) > 2) {
        err <- rep(Inf, kinase_num)
        for (i in 1:kinase_num) {
          k_temp <- k_set
          if (k_set[i] == 1) {
            err[i] <- 0
            k_temp[i] <- 0
            timma_temp <- timmaBinary(profile_data[, which(k_temp == 1)], sens, loo)
            err[i] <- sum(timma_temp$error)
            
          }
        }
        dummy <- min(err, na.rm = TRUE)
        ind_worst <- which.min(err)
        if (dummy < J[length(which(k_set == 1)) - 1]) {
          # if better result found !!J[step-1]!!
          k_set[ind_worst] <- 0
          if(verbosity){
            cat("Continuing excluding target #", ind_worst, "MAE = ", dummy, "\n")
          }
          
        } else {
          break
          if(verbosity){
            cat("break \n")
          }
          
        }
      }
      change <- 0
    }
  }
  
  temp <- which(J == min(J[which(J > 0)]))
  if (length(temp) > 1) {
    k_set <- X[, temp[2]]
  } else {
    k_set <- X[, temp[1]]
  }
  
  k_selected <- which(k_set == 1)
  profile_filtered <- unique(profile_data[, k_selected], MARGIN=2)
  timma <- timmaModel(profile_filtered, sens, loo)
  k_selected <- match(dimnames(profile_filtered)[[2]],dimnames(profile_data)[[2]])
  
  return(list(timma = timma, k_sel = k_selected))
} 
