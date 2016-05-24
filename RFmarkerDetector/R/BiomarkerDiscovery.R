


#' AUC multiple cross-validation
#' 
#' This function implements the AUCRF algorithm for identifying the variables (metabolites) most relevant
#' for the classification task 
#' @param data a n x p dataframe used to execute the AUCRF algorithm and perform a repetead CV of the AUCRF process.
#' The dependent variable must be a binary variable defined as a \code{factor} and codified as 0 for negatives (e.g controls)
#' and 1 for positivies (e.g. cases) 
#' @param seed a numeric value to set the seed of R's random number generator
#' @param ref_level the class assumed as reference for the binary classification
#' @param auc_rank the importance measure provided by \code{randomForest} for ranking the variables. There are
#' two options: MDG (default) and MDA 
#' @param auc_ntree the number of tree of each random forest model used 
#' @param auc_nfolds the number of folds in cross validation. By default a 5-fold cross validation is performed
#' @param auc_pdel the fraction of variables to be removed at each step. If \eqn{auc_pdel = 0}, it will be removed only one variable at each step
#' @param auc_colour the color chosen
#' @param auc_iterations a numeric that represents the number of cross validation repetitions    
#' @details Exploting the AUCRF algorithm, the fuction allows to identify the best performing 'parsimonious' model  
#' in terms of OOB-AUC and the most relevant variables (metabolites) involved in the prediction task. 
#' @import AUCRF
#' @importFrom stats relevel sd
#' @references Calle ML, Urrea V, Boulesteix A-L, Malats N (2011) 'AUC-RF: A new strategy for genomic pro-
#' filing with Random Forest'. Human Heredity
#' @export
#' @examples 
#' ## data(cachexiaData)
#' ## aucMCV(cachexiaData, ref_level = 'control')
#' @author Piergiorgio Palla 



aucMCV <- function(data, seed = 1234, ref_level = levels(data[, 2])[1], auc_rank = "MDG", auc_ntree = 500, auc_nfolds = 5, 
    auc_pdel = 0.2, auc_colour = "grey", auc_iterations = 5) {
    
    START <- proc.time()
    vip_measures <- c("MDA", "MDG")
    dataset = data
    names(dataset)[2] <- "class"
    
    
    set.seed(seed)
    
    if (!is.null(ref_level)) {
        ref_level = ref_level
        labels <- levels(data[, 2])
        if (!(ref_level %in% labels)) {
            stop("A problem occurred in opt: check ref_level parameter")
        }
        
        dataset[, 2] <- relevel(dataset[, 2], ref = ref_level)
    }
    
    levels(dataset[, 2]) = c(0, 1)
    ### NR:0; R:1 or C:0; R:1 or C:0; NR:1
    
    ## Here we perform RF-AUC feature selection on the entire data set: We use both MDG and MDA variable importance methods.
    
    # MDA
    if (!(auc_rank %in% vip_measures)) {
        auc_rank = "MDA"
    }
    
    
    rf.auc1 <- do.call("AUCRF", list(class ~ ., dataset[, 2:ncol(dataset)], pdel = auc_pdel, ranking = auc_rank, ntree = auc_ntree))
    print(rf.auc1)
    
    do.call("plot", list(rf.auc1, col = auc_colour))
    do.call("plot", list(rf.auc1, "ranking", col = "red"))
    
    cat(paste("... Performing a repeated ", auc_nfolds, "-folds CV\n ", sep = ""))
    cat(paste("Number of repetitions: ", auc_iterations, "\n", sep = ""))
    
    
    rf.aucv1 <- do.call("AUCRFcv", list(rf.auc1, nCV = auc_nfolds, M = auc_iterations))
    cat("...done!")
    
    ### there we plot the frequencies of the variable selected in the 5-fold CV
    
    plotVarFreq(varFrequency = rf.aucv1$Psel, thd = 0.55)
    
    ELAPSED <- proc.time() - START
    cat(paste("\nElapsed time:", round(ELAPSED[3], 2), "sec", sep = " "), "\n")
    invisible(rf.aucv1)
    
}


#' Variable Frequency Plot
#' 
#' This function plots the probability of selection of each variable (metabolite) as the proportion of times that is selected
#' by AUCRF method
#' @param varFrequency a numeric vector with the probabilities of selection of each input variable 
#' @param thd a numeric value that indicates the lower limit of probability of variables to be represented. The default
#' value is 0.4 (40\%)
#' @param color the color of the graph
#' @param hcol the color of the horizontal line, indicating the lower limit of the probability of selection of variables
#' @export
#' @examples
#' ## data(cachexiaData)
#' ## rf.aucv1 <- aucMCV(cachexiaData, ref_level = 'control')
#' ## plotVarFreq(varFrequency = rf.aucv1$Psel,  thd=0.55)
#' @author Piergiorgio Palla
#' 
plotVarFreq <- function(varFrequency, thd = 0.4, color = "blue", hcol = "red") {
    
    ##### Usually our vector will be a named vector in which the name of each element will be the var name To allow a comparision
    ##### between different approaches we sort the names and thus the vector in alphabetical order
    
    if (!is.null(names(varFrequency))) {
        unordered_names <- names(varFrequency)
        ordered_names <- sort(unordered_names)
        varFrequency <- varFrequency[ordered_names]
    }
    
    topFrequenciesIdx <- which(varFrequency > thd)
    cols <- rep("gray", length(varFrequency))
    xarguments <- rep("", length(varFrequency))
    
    cols[topFrequenciesIdx] = color
    xarguments[topFrequenciesIdx] = names(topFrequenciesIdx)
    xarguments_full_names <- xarguments
    xarguments <- abbreviate(xarguments, minlength = 5)
    cat(length(xarguments))
    cat(xarguments)
    
    par(mar = c(5, 5, 2, 2))
    
    b = barplot(varFrequency * 100, col = cols, ylim = c(0, 105), names.arg = xarguments, ylab = "Candidate biomarker frequencies (%)", 
        space = 0.7, cex.names = 0.75)
    idx = 1:length(varFrequency)
    lab = rep("", length(varFrequency))
    
    text(x = b, y = varFrequency * 100 + 4, labels = xarguments)
    abline(h = thd * 100, col = hcol, lty = 3)
    
    
    
} 
