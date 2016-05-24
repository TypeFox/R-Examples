#' Filter targets
#' 
#' A function to filter targets based on their corration with the drug sensitivity
#' @param profile drug-target interaction data
#' @param sens drug sensitivity data
#' @param sp an integer to specify the starting point for the sffs search algorithm. The number cannot be larger than the total number of targets in the drug-target interaction data. 
#' By default, the starting point is the first target, namely, sp = 1.
#' @param max_k an integer to specify the maximal number of targets that can be selected by the sffs algorithm. In practice it is advised to keep it under 10 
#' as the number of sensitivities to be predicted will increase exponentially. By default, max_k = 2.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
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
#' result<-floating2(tyner_interaction_binary, tyner_sensitivity[,1], sp = 1, max_k = 5)
#' }
floating2 <- function(profile, sens, sp = 1, max_k = 2, verbosity = FALSE) {
    drug_num <- nrow(profile)
    k_num <- ncol(profile)
    loo <- 1
    corrList <- cor(profile, sens)
    keep_list <- which(corrList > 0)
    new_profile <- profile[, keep_list]
    k <- length(keep_list)
    # get the new_initial_list
    sp_new <- which(keep_list == sp)
    if (length(sp_new) != 0) {
        new_initial_list <- rep(0, k)
        new_initial_list[sp_new] <- 1
    } else {
        new_initial_list <- c(1, rep(0, k - 1))
        cat("start from 1 not sp")
    }
    return(sffsBinary2(new_profile, sens, sp_new, max_k, loo, new_initial_list, verbosity))
} 
