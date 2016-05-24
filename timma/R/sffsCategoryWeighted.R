#' Model selection with sffs for the multi-class drug-target interaction data using one.sided and weighted TIMMA model
#' 
#' A function to select the most predictive targets with sffs for the multi-class drug-target interaction data using
#' one.sided and weighted TIMMA model
#' 
#' @param profile_data drug-target interaction data which is a matrix with drugs as row indexes and targets
#' as column indexes.
#' @param sens a drug sensitivity vector.
#' @param sp an integer to specify the starting point for sequential forward floating search (sffs) search 
#' algorithm to navigate the target set space. By default, sp = 1.
#' @param max_k an integer to sepcify the maximum number of targets that can be selected by the sffs
#' algorithm. By default, max_k = 2. In practice it should not be over than 10 as the number of target combinations will increase exponentially.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#' @param class  an integer number to specify the number of classes in the drug-target interaction data. 
#' For a binary drug-target interaction data, class = 2. For a multi-class drug-target interaction data, 
#' class should be the number of classes. 
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' @return A list containing the following components:
#' \item{timma}{a list contains: the predicted efficacy matrix, prediction error and predicted drug sensitivity}
#' \item{k_sel}{the indexes for selected targets}
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_multiclass)
#' data(tyner_sensitivity)
#' results<-sffsCategoryWeighted(tyner_interaction_multiclass, tyner_sensitivity[, 1], class = 6)
#' }

sffsCategoryWeighted <- function(profile_data, sens, sp = 1, max_k = 2, loo = TRUE, class, verbosity = FALSE) {
    
    # get the number of the drugs
    drug_num <- nrow(profile_data)
    # get the number of the kinases
    kinase_num <- ncol(profile_data)
    
    # binary vector to indicate if a kinase is in the cancer-specific target sets
    new_initial_list <- rep(0, kinase_num)
    # set the sp kinase indicator to be 1
    new_initial_list[sp] <- 1
    
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
            # timma_temp<-TIMMA(profile_temp, sens, loo, step, drug_num)
            timma_temp <- timmaCategoryWeighted(profile_temp, sens, loo, class)
            # the score for X_k, the lower score is better
            J[step] <- mean(timma_temp$error)
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
            
            # cat(sprintf('k=%d, error=%4.2f', step, J[step]))
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
            
            # step 1: Inclusion space<-search_space1(drug_num, k_set, profile_data, sens)
            for (i in 1:kinase_num) {
                if (k_set[i] == 0) {
                  # including new kinase i
                  err[i] <- 0
                  set_add_one <- k_set
                  set_add_one[i] <- 1
                  err1 <- timmaCategoryWeighted(profile_data[, which(set_add_one == 1)], sens, loo, class)$error
                  # err1<-TIMMA1(profile_data[,c(which(k_set==1), i)], sens, loo, sum(k_set)+1, drug_num)$error
                  err[i] <- mean(err1)
                }
            }
            
            # select X_{k+1} the most significant feature with respect to X_k
            dummy <- min(err, na.rm = TRUE)
            ind1 <- which.min(err)
            k_set[ind1] <- 1
            # sprintf('Inclusion kinase=%d error=%4.2f\n', ind1, dummy)
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
                  # timma_temp<-TIMMA(profile_temp, sens, loo, sum(k_temp), drug_num)
                  timma_temp <- timmaCategoryWeighted(profile_temp, sens, loo, class)
                  err[i] <- mean(timma_temp$error)
                }
            }
            worst_err <- min(err)
            ind_worst <- which.min(err)
            if (err[ind1] != worst_err && worst_err < J[step]) {
                # the new added kinase is not the least significant feature
                k_set[ind_worst] <- 0  # exclude the least significant feature
                # cat(sprintf('Exclusion kinase=%d error=%4.2f/n', ind_worst, worst_err))
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
                    # timma_temp<-TIMMA(profile_data[,which(k_temp==1)], sens, loo, sum(k_temp), drug_num)
                    timma_temp <- timmaCategoryWeighted(profile_data[, which(k_temp == 1)], sens, loo, class)
                    err[i] <- mean(timma_temp$error)
                    
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
    timma <- timmaCategoryWeighted(profile_filtered, sens, loo, class)
    k_selected <- match(dimnames(profile_filtered)[[2]],dimnames(profile_data)[[2]])
    return(list(timma = timma, k_sel = k_selected))
} 
