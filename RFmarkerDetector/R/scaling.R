
#' Pareto scaling method performed on the columns of the data table (i.e. metabolite concentrations measured by 1H NMR
#' or binned 1H NMR spectra)
#' 
#' The function provides a data pretreatment approach called Pareto Scaling. Each column of the table is given a mean of zero by substracting the column
#' column mean from each value in the column; then each value in each column is divided by a scaling factor, represented by the square 
#' root of the standard deviation of the column values.
#' 
#' @param data a n x p matrix of n observations and p predictors. If the first two columns of the matrix
#' represent respectively the sample names and the class labels associated to each sample, the scaling method should not
#' include these two columns
#' @param exclude a boolean variable. If set to True the scaling method will exclude the first two columns. 
#' @return a scaled version of the input matrix
#' @author Piegiorgio Palla
#' @details 
#' This function is useful when variables have significantly different scales. It is  generally the preferred option 
#' in NMR Metabolomics because it is a good compromise between no scaling (centering) and auto scaling
#' @examples
#' #' ## load the included example dataset
#' data(cachexiaData)
#' ## call paretoscale with the parameter exclude set to TRUE (default) 
#' ## in order to exclude the first two columns of the dataset from scaling
#' data.scaled <- paretoscale(cachexiaData, exclude = TRUE) 
#' @export

paretoscale <- function(data, exclude = T) {
    
    if (exclude == T) {
        # Here we extract numeric data and perform Pareto scaling
        sample_classes <- data[, 1:2]
        x <- data[, 3:dim(data)[2]]
    } else {
        sample_classes <- NULL
        x <- data
    }
    # Here we perform centering
    x.centered <- apply(x, 2, function(x) x - mean(x))
    # Then we perform scaling on the mean-centered matrix
    x.sc <- apply(x.centered, 2, function(x) x/sqrt(sd(x)))
    x.sc <- cbind(sample_classes, x.sc)
    
}

#' Unit variance scaling method performed on the columns of the data (i.e. metabolite concentrations measured by 1H NMR
#' or binned 1H NMR spectra)
#' 
#' The function provides a data pretreatment approach called Autoscaling (also known as unit variance scaling). The data for each variable (metabolite)
#' is mean centered and then divided by the standard deviation of the variable. This way each variable will have zero mean and unit standard deviation.  
#' 
#' @param data a n x p data frame with n observations and p columns. While the first two columns usually represent the names of the samples and the 
#' class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra
#' @param exclude a boolean variable which stores a simple True/ False setting. If set to True the scaling method will exclude the first two columns. 
#' @return a scaled version of the input matrix
#' @author Piegiorgio Palla
#' @examples
#' ## load the included example dataset
#' data(cachexiaData)
#' ## call autoscale with the parameter exclude set to TRUE (default) 
#' ## in order to exclude the first two columns of the dataset from scaling
#' data.scaled <- autoscale(cachexiaData, exclude = TRUE)
#' @export

autoscale <- function(data, exclude = T) {
    
    if (exclude == T) {
        # Here we extract numeric data and perform Pareto scaling
        sample_classes <- data[, 1:2]
        x <- data[, 3:dim(data)[2]]
    } else {
        sample_classes <- NULL
        x <- data
    }
    # Here we perform autoscaling, centering each column and dividing it by its standard deviation
    x.sc <- scale(x, center = T, scale = T)
    x.sc <- cbind(sample_classes, x.sc)
    x.sc
    
}

#' Mean centering performed on the columns of the data (i.e. metabolite concentrations measured by 1H NMR
#' or binned 1H NMR spectra)
#' 
#' The function allows to have each predictor (column) centered on zero. The average value of each predictor is substracted to each value in the column.
#' 
#' @param data a n x p data frame with n observations and p columns. While the first two columns usually represent the names of the samples and the 
#' class labels related to each sample respectively, the remaining columns represent metabolite concentrations measured by 1H NMR or bins of 1H NMR spectra
#' @param exclude a logical variable which stores a simple True / False setting. If set to True the scaling method will exclude the first two columns. 
#' @return a mean centered version of the input matrix
#' @author Piergiorgio Palla
#' @examples
#' ## load the included example dataset
#' data(cachexiaData)
#' ## call meanCenter with the parameter exclude set to TRUE (default) 
#' ## in order to exclude the first two columns of the dataset from scaling
#' data.scaled <- meanCenter(cachexiaData, exclude = TRUE)
#' @export
meanCenter <- function(data, exclude = T) {
    
    if (exclude == T) {
        # Here we extract numeric data and perform Pareto scaling
        sample_classes <- data[, 1:2]
        x <- data[, 3:dim(data)[2]]
    } else {
        sample_classes <- NULL
        x <- data
    }
    # Here we perform mean centering
    x.sc <- scale(x, center = T, scale = F)
    x.sc <- cbind(sample_classes, x.sc)
} 
