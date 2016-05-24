#' SFFS switch function
#' 
#' A function to choose which sffs function to run. There are six sffs algorithms for choosing. 
#' 
#' @param profile_data drug-target interaction data which is a matrix with drugs as row indexes and targets
#' as column indexes.
#' @param sens a drug sensitivity vector.
#' @param sp an integer to specify the starting point for the sffs search algorithm. The number cannot 
#' exceed the total number of targets in the drug-target interaction data. By default, the starting 
#' point is the first target, namely, sp = 1.
#' @param max_k an integer to sepcify the maximum number of targets that can be selected by the sffs
#' algorithm. By default, max_k = 2. In practice it should not be over than 10 as the number of target combinations will increase exponentially.
#' @param loo a logical value indicating whether to use the leave-one-out cross-validation in the model
#' selection process. By default, loo = TRUE. 
#' @param class  an integer to specify the number of classes in the drug-target interaction data. 
#' For a binary drug-target interaction data, class = 2. For a multi-class drug-target interaction
#'  data, class should be the number of classes. 
#' @param averaging a parameter to specify which one of the averaging algorithms will be applied 
#' in the model construction. By default, averaging = "one.sided", which is the original model 
#' construction algorithm. When averaging = "two.sided", a modified averaging algorithm will be used.
#' These two variants only differ for the case where the minimization and maximization rules are not
#' simultaneously satisfied. For example, for a queried target set if the supersets but not the subsets
#' can be found in the training data, the one.sided algorithm will take the prediction from the averages
#' on the supersets sensitivities using the minimization rule. The two.sided algorithm, however, will 
#' lower the predicted sensitivity by averaging it with 0, which is the theoretical lower boundary of 
#' the sensitivities that could be obtained in the subsets. 
#' @param weighted a parameter to specify if the similarity between the queried target set and 
#' its subsets/supersets is considered as a weight factor in the averaging. When weighted =T RUE,
#' the similarity is considered as a weight factor such that those related target sets will be
#' weighted more in the final predictions.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' 
#' @return A list containing the following components:
#' \item{timma}{a list contains: the predicted efficacy for target combinations, prediction error and predicted drug sensitivity}
#' \item{k_sel}{the indexes for selected targets}
#' @author Liye He \email{liye.he@@helsinki.fi} 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' results<-sffs(tyner_interaction_binary, tyner_sensitivity[, 1], max_k = 8)
#' }


sffs <- function(profile_data, sens, sp = 1, max_k = 2, loo = TRUE, class = 2, averaging = "one.sided", weighted = FALSE, verbosity = FALSE) {
    # nbits: if it is 2, it calls the binary timma; otherwise, it calls the category timma. 

    if (averaging == "one.sided") {
        if(class == 2){
          return(sffsBinary(profile_data, sens, sp, max_k, loo, verbosity))
        }
        else if(class > 2 && weighted == FALSE){
          return(sffsCategory(profile_data, sens, sp, max_k, loo, class, verbosity))
        }else{
          return(sffsCategoryWeighted(profile_data, sens, sp, max_k, loo, class, verbosity))
        }
    } else if (averaging == "two.sided") {
      if(class == 2){
        return(sffsBinary1(profile_data, sens, sp, max_k, loo, verbosity))
      }
      else if(class > 2 && weighted == FALSE ){
        return(sffsCategory1(profile_data, sens, sp, max_k, loo, class, verbosity))
      }else{
        return(sffsCategoryWeighted1(profile_data, sens, sp, max_k, loo, class, verbosity))
      }
    } 
} 
