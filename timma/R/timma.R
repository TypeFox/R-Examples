#' Main function for the timma package
#' 
#' Target inhibition inference using maximization and minimization averaging
#'
#' @param x a drug-target interaction matrix. Row names are drug names and column names are target names.
#' @param y a normalized drug sensitivity vector.
#' @param sp an integer to specify the starting point for the sffs search algorithm. The number cannot be larger than the total number of targets in the drug-target interaction data. 
#' By default, the starting point is the first target, namely, sp = 1.
#' @param max_k an integer to specify the maximal number of targets that can be selected by the sffs algorithm. In practice it is advised to keep it under 10 
#' as the number of sensitivities to be predicted will increase exponentially. By default, max_k = 5.
#' @param filtering a logical parameter to determine whether the targets should be filtered before the model selection. 
#' By default, the value is FALSE, meaning that all the available targets will be considered in the model selection. 
#' If the value is TRUE, those targets that are negatively correlated with the drug sensitivity data will be removed.
#' @param class an integer to specify the number of classes in the drug-target interaction data. For a binary drug-target 
#' interaction data, class = 2. For a multi-class drug-target interaction data, class should be the number of classes.
#' @param averaging a parameter to specify which one of the averaging algorithms will be applied in the model construction. 
#' By default, averaging = "one.sided", which is the original model construction algorithm. When averaging = "two.sided", 
#' a modified averaging algorithm will be used. These two variants only differ for the case where the minimization and 
#' maximization rules are not simultaneously satisfied. For example, for a queried target set if the supersets but not the 
#' subsets can be found in the training data, the one.sided algorithm will take the prediction from the averages on the 
#' supersets sensitivities using the minimization rule. The two.sided algorithm, however, will lower the predicted 
#' sensitivity by averaging it with 0, which is the theoretical lower boundary of the sensitivities that could be 
#' obtained in the subsets.
#' @param weighted When averaging = "weighted", the similarity between the queried target set and 
#' its subsets/supersets is considered as a weight factor in the averaging, such that those related target 
#' sets will be more weighted in the final predictions.
#' @param verbosity a boolean value to decide if the information should be displayed. If it is TRUE, the information
#' will be displayed while the model is running. Otherwise, the information will not be displayed. By default, it is
#' FALSE.
#' @param use When use = "observed", the true drug sensitivity data will be used for drawing target inhibition network.
#' When use = "predicted", the predicted drug sensitivity data will be used for drawing target inhibition network.
#' @return an R image of the input and output data. 
#' @author Jing Tang \email{jing.tang@@helsinki.fi} 
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' 
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' median_sensitivity<-tyner_sensitivity[, 1]
#' results<-timma(tyner_interaction_binary, median_sensitivity)
#' }

timma <- function(x, y, sp = 1, max_k = 5, filtering = FALSE, class = 2, averaging = "one.sided", weighted = FALSE, verbosity = FALSE, use = "observed") {
  
  # get drug names
  drug_names <- dimnames(x)[[1]]
  
  # get kinase names
  kinase_names <- dimnames(x)[[2]]
  
  cat("----------------Start Running TIMMA----------------------------- \n")
  if (filtering == FALSE) {
    # sffs
    float<-sffs(x, y, sp, max_k, loo = TRUE, class, averaging, weighted, verbosity)
  } else {
    float <- floating2(x, y, sp, max_k, verbosity)
  }
  cat("----------------Complete running TIMMA model-------------------- \n")
  
  err <- float$timma$error
  error <- mean(err)
  RMSE <- sqrt(mean(err^2))
  RMSE_baseline <- sqrt(mean((y - median(y))^2))
  
  # write selected Kinase file
  x <- data.frame(x)
  k_select <- float$k_sel
  select_kinase_names <- findSameSet(x, k_select, kinase_names)
  #select_file_name <- paste("selectedKinase", model, ".csv", sep = "")
  select_file_name <- "selectedTargets.csv"
  profile_select <- x[, k_select]
  
  if (max_k > 1) {
    #dimnames(profile_select)[[2]] <- select_kinase_names  
    select_target <- cbind(profile_select, y, float$timma$prediction)
    dimnames(profile_select)[[2]] <- select_kinase_names
    dimnames(select_target)[[2]] <- c(select_kinase_names, "sensitivity", "LOO sensitivity")
    write.table(select_target, file = select_file_name, sep = ",", col.names = NA)
  } else {
    profile_select <- rbind(select_kinase_names, profile_select)
  }
  cat("----------------Saving selectedTargets.csv---------------------- \n")
  
  # write Timma
  if (class == 2) {
    gc_timma <- graycode3(length(k_select))
    gc_names <- graycodeNames(length(k_select), select_kinase_names, gc_timma$gc_row, gc_timma$gc_col)
    #timma_file <- paste("Timma", model, ".csv", sep = "")
    timma_file <- "predictedSensitivities.csv"
    nr <- gc_names$nr
    nc <- t(gc_names$nc)
    timma_row <- nrow(nr) + nrow(nc)
    timma_col <- ncol(nr) + ncol(nc)
    timma <- array("", dim = c(timma_row, timma_col))
    timma[(nrow(nc) + 1):timma_row, 1:ncol(nr)] <- nr
    timma[1:nrow(nc), (ncol(nr) + 1):timma_col] <- nc
    timma[(nrow(nc) + 1):timma_row, (ncol(nr) + 1):timma_col] <- float$timma$dummy
    write.table(timma, file = timma_file, sep = ",", col.names = FALSE, row.names = FALSE)
    cat("----------------Saving predictedSensitivities.csv--------------- \n")
    
    target_comb_rank <- targetRank(profile_select, timma)
    write.table(target_comb_rank, file = "predictedTargetScoring.csv",sep=",",row.names = FALSE)
    cat("----------------Saving predictedTargetScoring.csv--------------------- \n")
    
    drug_comb_rank <- drugRank(profile_select, timma, y)
    write.table(drug_comb_rank, file = "predictedDrugScoring.csv", sep = ",", row.names = FALSE)
    cat("----------------Saving predictedDrugScoring.csv--------------------- \n")
    
    # write the R image
    save(x,profile_select, timma, y, file="result.RData")
    cat("----------------Saving result.RData---------------------- \n")
    
    #loo_prediction <- float$timma$prediction
    if(use == "observed"){
      loo_prediction <- y
    }else if(use == "predicted"){
      loo_prediction <- float$timma$prediction
    }else{
      stop("use must be either observed or predicted")
    }
    
    one<-which(loo_prediction>0.5)
    zero<-which(loo_prediction<=0.5)
    SENS<-loo_prediction
    SENS[one]<-1
    SENS[zero]<-0
   
    draw_data<-cbind(profile_select, SENS)
    drawGraph(draw_data)
    cat("----------------Saving targetInhibitionNetwork.pdf-------------- \n")
    cat("----------------Saving targetInhibitionNetwork.nnf-------------- \n")
  }
  
  current_dir <- getwd()
  cat("Analysis finished. All the results are saved in", current_dir)

  return(list(profile_input = x, profile_select = profile_select, prediction = timma, sensitivity = y))
} 
