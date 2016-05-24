#' Filtering 'low quality' variables from the original dataset 
#' 
#' This function takes the original dataset and filters those variables with a definite (usually relevant) percentage of zero-values
#' @param data a n x p data frame with n observations and p columns. While the first two columns usually represent the names of the samples and the 
#' class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra
#' @param exclude a logical variable which stores a simple True / False setting. If set to True the filtering method will exclude the first two columns.
#' @param threshold the percentage of zero values of a variable above which it will be eliminated from the dataset (default: 0.50)
#' @return a list containing the filtered dataset, a vector with the names of the varables excluded and a vector with the indexes of the variables eliminated
#' @examples 
#'  ## load the included example dataset
#' data(cachexiaData)
#' ## call lqvarFilter with the parameter exclude set to TRUE (default) 
#' ## in order to exclude the first two columns of the dataset from scaling
#' res <- lqvarFilter(cachexiaData, threshold = 0.4, exclude = TRUE)
#' data.filtered <- res$filtered_dataset 
#' @author Piergiorgio Palla
#' @export

lqvarFilter <- function(data, threshold = 0.5, exclude = T) {
    
    if (exclude == T) {
        sample_classes = data[, 1:2]
        input <- data[, 3:ncol(data)]
    } else {
        sample_classes = NULL
        input <- data
    }
    
    n_of_els <- dim(input)[1]
    
    idx = c()
    for (j in 1:ncol(input)) {
        
        current_col = input[, j]
        negative_idx <- which(current_col <= 0)
        n_of_negatives <- length(negative_idx)
        
        ratio <- n_of_negatives/n_of_els
        if (ratio >= threshold) {
            idx <- append(idx, j)
        }
        
    }
    
    if (!(is.null(idx) || length(idx) == 0)) {
        excluded_vars <- names(input)[idx]
        filtered_input <- input[, -idx]
    } else {
        excluded_vars = NULL
        filtered_input = input
    }
    
    
    reduced_dataset <- cbind(sample_classes, filtered_input)
    
    return(list(filtered_dataset = reduced_dataset, excluded_vars = excluded_vars, idx = idx))
    
}

#' Computing relative standard deviation of a vector
#' 
#' This function computes the relative standard deviation (also known as coefficient of variation) of a numeric vector defined
#' as the ratio of the standard deviation to the mean of the vector elements, expressed as percentage
#' @param v a numeric vector
#' @return the value of the coefficient of variation of the input vector expressed as a percentage and rounded to two
#' decimal places
#' @author Piergiorgio Palla
#' @details the coefficient of variation shows the extent of variability in relation to mean of the population.
#' It is expressed as a percentage. Lower values indicate lower variability.
#' @examples
#' v <-  runif(10, min = 5, max = 30)
#' rsd(v) 
#' @export

rsd <- function(v) {
    
    std <- sd(v)  #standard deviation
    mean_value <- mean(v)  # mean
    
    
    rsd <- (std/mean_value) * 100  # relative standard deviation expressed as a percentage
    return(round(rsd, digits = 2))
    
}


#' Filtering less informative variables
#' 
#' \code{rsdFilter} removes from the dataframe the predictors with a relative standard deviation less
#' than or equal to an inserted threshold
#' @param data a n x p data frame with n observations and p columns. While the first two columns usually represent the names of the samples and the 
#' class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra
#' @param threshold a numeric value representing a limit: each predictor with a relative standard deviation lower than that 
#' will be removed form the dataframe
#' @param exclude a logical variable which stores a simple True / False setting. If set to True the filtering method will exclude the first two columns.
#' @return a list containing the filtered dataset, a vector with the names of the varables excluded and a vector with the indexes of the variables eliminated
#' @author Piergiorgio Palla
#' @examples
#'  ## load the included example dataset
#' data(cachexiaData)
#' ## call rsdFilter with the parameter exclude set to TRUE (default) 
#' ## in order to exclude the first two columns of the dataset from scaling
#' data.filtered <- rsdFilter(cachexiaData, threshold = 15, exclude = TRUE)
#' @export  

rsdFilter <- function(data, threshold, exclude = T) {
    
    if (exclude == T) {
        sample_classes = data[, 1:2]
        input <- data[, 3:ncol(data)]
    } else {
        sample_classes = NULL
        input <- data
    }
    
    ### Here we compute the relative std of each column We want to eliminate the 'invariant' variables
    
    rsd_vector <- apply(input, 2, rsd)
    
    ## Here we select the informative variables
    col_idx <- which(rsd_vector >= threshold)
    
    excluded_idx <- which(rsd_vector < threshold)
    
    filtered_input <- input[, col_idx]
    filtered_dataset <- cbind(sample_classes, filtered_input)
    excluded_vars = names(input)[-col_idx]
    
    return(list(filtered_dataset = filtered_dataset, excluded_vars = excluded_vars, idx = excluded_idx))
    
    
} 
