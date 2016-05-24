

#' Combinatorial Monte Carlo CV
#' 
#' This function performs a Monte Carlo CV for each of the Random Forest model grown considering
#' all the k-combinations of the n input variables of the original dataset, with k ranging from 2 to n.
#' It allows to get the most performing Random Forest model in terms of the AUC of the ROC curve and
#' to obtain the most relevant input variables (metabolites or bins) associated with it.
#' @param dataset a n x p dataframe used to build the models. The first two columns
#' must represent respectively the sample names and the class labels related to each sample
#' @param parameters a list including the following parameters: \itemize{
#' \item ntree the number of trees of each Random Forest model
#' \item nsplits the number of random splittings of the original dataset into training and test data sets
#' \item test_prop the percentage (expressed as a real number) of the observations of the original dataset
#' \item kmax the maximum number of inputs to combine. 
#' }
#' @details The function computes all the k-combinations of the n input variables, with k ranging from 2 to n.
#' Each combination corresponds to a dataset on which the function will grow a Random Forest model,
#' performing a Monte Carlo CV. Then it will provide the best performing model in terms of the AUC of the ROC curve
#' and the most relevant variables associated with it.
#' @importFrom utils combn
#' @return a list containing the most performing Random Forest model
#' #' @examples
#' ## data(cachexiaData) 
#' ## params <- list(ntrees = 100, nsplits = 10, test_prop = 1/3)
#' ## res <- combinatorialRFMCCV(dataset = cachexiaData[,1:10], parameters = params)
#' ## This task may take a long time depending on the 
#' ## dimension of the dataset and on the parameters provided
#' @export
#' @author Piergiorgio Palla 


combinatorialRFMCCV <- function(dataset, parameters = list(ntrees = 500, nsplits = 100, test_prop = 1/3, kmax = 5)) {
    
    START <- proc.time()
    
    original_data = dataset
    start_idx <- 3
    last_idx <- dim(original_data)[2]
    
    ## These are the column indexes of original data - i.e. concentration matrix {3,4,...13}
    col_idx <- seq(start_idx, last_idx)
    
    
    ### Test Params #### The label_order param states that 'cases' are 'positives' and 'ctrl' are negatives test_params <-
    ### list(ntree= 800, nsplits = 100, test_prop = 1/3, label_order=c('ctrl', 'case'))
    test_params <- parameters
    
    # combn(v,2)
    num_of_models = 1
    
    if (hasArg(parameters) && (!is.null(parameters$kmax)) && (kmax <= length(col_idx))) {
        kmax <- parameters$kmax
    } else {
        kmax = length(col_idx)
    }
    
    ## Here we consider the j-combinations from the col_idx array Each combination represents a set of metabolites. We will
    ## consider the datasets extracted form the original dataset considering the column indexes correpsonding to these
    ## combinations.  On each of these datasets, we will generate a RF model.  We are trying to find the best combination of
    ## metabolites. The best set of candidate biomarkers will be the one with the RF model with the highest performance on the
    ## corresponding dataset.  To compare each model we will use the Area under the Roc Curve (AUC) built considering the
    ## erformance of each RF model.
    top_models = list()
    top_metabolites = list()
    auc_values = c()
    cat("\n\n ... Performing CV\n\n ")
    ### We are collecting the k-combinations from the inputa data index, with k ranging from 10 to 2
    for (j in 2:(kmax - 1)) {
        cat(paste(j, " - combinations\n"))
        combinations <- combn(col_idx, j)  # a j x n_of_combinations matrix
        num_of_models = num_of_models + dim(combinations)[2]
        ### Here we try to find the best p
        
        res <- getBestRFModel(combinations, dataset, test_params)
        auc_values <- append(auc_values, res$auc)
        top_models[[j - 1]] <- res$best_model_set
        top_metabolites[[j - 1]] <- res$biomarker_set
        
        ####### Accuracy Recall x Model #####
        
        # perf <- rfMCCVPerf(top_models[[j-1]]$mod)
        
        # cat(paste('Accuracy: ', perf$avg_accuracy)) cat('\n') cat(paste('Recall: ', perf$avg_recall)) cat('\n')
        
        
        
    }
    
    ### Best Performing Model ###
    max_idx <- which.max(auc_values)
    print(paste("Max AUC: ", max(auc_values)))
    top_model <- top_models[[max_idx]]
    top_metabolite_set <- top_metabolites[[max_idx]]
    
    cat("\n Best Performing model:\n")
    # print(top_model$mod)
    
    cat("\nBest combination of Metabolites: ")
    cat(names(original_data)[top_metabolite_set])
    cat("\n")
    cat(paste("size", length(top_metabolite_set)))
    cat("\n")
    
    
    ## Average Accuracy and Recall
    
    
    perf <- rfMCCVPerf(top_model$mod)
    
    cat(paste("Accuracy: ", perf$avg_accuracy))
    cat("\n")
    cat(paste("Recall: ", perf$avg_recall))
    cat("\n")
    
    
    ### Plot averaged ROC curve and AUC confidence interval
    plot.mccv(top_model)
    
    
    
    
    
    ELAPSED <- proc.time() - START
    cat(paste("Elapsed time:", round(ELAPSED[3], 2), "sec", sep = " "), "\n")
    
    return(list(best_model = top_model))
    
} 
