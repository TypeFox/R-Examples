




#===============================================================================
# Filename: checkDataQuality.R
#
# Purpose: Check data quality and generate summary statistics
# Authors: Madhav Kumar & Shreyes Upadhyay
# Date: 15Sep2013
# Version: 1.0
# Dependencies: None
# Packages required: None
#
# Inupts: 1. R Data Frame; no default
#         2. Csv output filename for numeric variables
#         3. Csv output filename for categorical variables
#         4. Numeric cutoff - Minimum number of unique values needed for a 
#            to be considered numeric; default = -1 (See notes for details)
#
# Outputs: 1. Numeric: Csv file containing data quality report (number missing, 
#             number unqiue etc. and summary statistics)
#          2. Categorical: Csv file containing data quality report and 
#             counts of top 10 most frequent categories for each variable
#
# Limitations: 1. Only works for numeric, character, and factor variable 
#                 types in R. Does not recognize date and boolean
#===============================================================================

#' @export checkDataQuality

#==========================================
# Data quality checks
#==========================================
checkDataQuality <- function(data,  
                             out.file.num, 
                             out.file.cat,
                             numeric.cutoff= -1){
  options(scipen = 999)
  start.time <- Sys.time()
  
  # columns in data set
  cols <- 1:ncol(data)
  
  # determine which variables are categorical, i.e., either character or numeric
  cats <- sapply(cols,  function(i) is.factor(data[, i]) || is.character(data[,i]) || 
                   length(unique(data[, i])) <= numeric.cutoff)
  cats <- which(cats == TRUE)
  
  # determine which variables are numeric
  nums <- sapply(cols,  function(i) is.numeric(data[, i]) & 
                   length(unique(data[, i])) > numeric.cutoff)
  nums <- which(nums == TRUE)
  
  # supplemental functions
  maxNA <- function(x){
    if (all(is.na(x))){
      return(NA)
    } else return(max(x, na.rm= TRUE))
  }
  
  minNA <- function(x){
    if (all(is.na(x))){
      return(NA)
    } else return(min(x, na.rm= TRUE))
  }
  
  #========================================
  # Quality check on numeric variables
  #========================================
  
  if (length(nums) > 0){
    # take data with only numeric variables
    num.data <- data[, nums]
    
    # non-missing values
    n.non.miss <- colSums(!is.na(num.data))
    
    # missing values
    n.miss <- colSums(is.na(num.data))
    
    # missing percentage
    n.miss.percent <- 100*n.miss/nrow(num.data)
    
    # unique values
    n.unique <- apply(num.data, 2, unique)
    
    # count of unique
    n.unique <- simplify2array(lapply(n.unique, length))  
    
    # mean value
    n.mean <- apply(num.data, 2, mean, na.rm= TRUE)
    
    # minimum value
    n.min <- apply(num.data, 2, minNA)
    
    # maximum value
    n.max <- apply(num.data, 2, maxNA)
    
    # quantiles
    n.quant <- apply(num.data, 2, quantile, probs= c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 
                                                     0.9, 0.95, 0.99), na.rm= TRUE)
    
    # combine all results in data frame
    n.output <- rbind(n.non.miss, n.miss, n.miss.percent, n.unique, 
                      n.mean, n.min, n.quant, n.max)
    
    # transpose output data 
    n.output <- data.frame(t(n.output))
    
    # round results to two decimal digits
    n.output <- round(n.output, 2)
    
    # add col names to output data
    names(n.output) <- c('non-missing', 'missing', 'missing percent', 
                         'unique', 'mean', 'min', 'p1', 'p5',  'p10', 
                         'p25', 'p50', 'p75', 'p90', 'p95', 'p99', 'max')
    write.csv(n.output, out.file.num, row.names= TRUE)
    cat('Check for numeric variables completed')
    cat(' // ')
    cat('Results saved to disk') 
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(' // ')
  print(time.taken)
  
  
  #========================================
  # Quality check on categorical variables
  #========================================
  start.time <- Sys.time()
  
  if (length(cats) > 0){
    # take data with only categorical variables
    cat.data <- data[, cats]
    
    # convert all to character
    cat.data[, 1:ncol(cat.data)] <- lapply(cat.data[, 1:ncol(cat.data)], as.character)
    
    # set all missing values to NA
    cat.data[, 1:ncol(cat.data)] <- lapply(cat.data[, 1:ncol(cat.data)],
                                           function(x){
                                             ifelse(x == "", NA, x)
                                           })
    
    # non-missing values
    n.non.miss <- colSums(!is.na(cat.data))
    
    # missing values
    n.miss <- colSums(is.na(cat.data))
    
    # missing percentage
    n.miss.percent <- 100*n.miss/nrow(cat.data)
    
    # unique values
    n.unique <- apply(cat.data, 2, unique)
    
    # count of unique
    n.unique <- simplify2array(lapply(n.unique, length))
    
    # combine all results in data frame
    n.output <- rbind(n.non.miss, n.miss, n.miss.percent, n.unique)
    
    # transpose output data 
    n.output <- data.frame(t(n.output))
    
    # round results to two decimal digits
    n.output <- round(n.output, 2)
    
    # categories and their frequencies
    n.categories <- apply(cat.data, 2, function(x) sort(table(x), decreasing= TRUE))
    
    # count maximum categories
    max.cat <- max(unlist(lapply(n.categories, length)))
    if (max.cat > 10) max.cat <- 10
    
    # create empty variables in outupt
    cat.names <- paste(rep(c("cat", "freq"), max.cat), rep(1:max.cat, each= 2), sep= "_")
    n.output[, cat.names] <- ""
    
    n.output <- lapply(row.names(n.output), 
                       function(x){
                         tmp <- n.output[row.names(n.output) == x,]
                         freqs <- length(n.categories[[x]])
                         freqs <- pmin(10, freqs)
                         if (length(freqs) == 1 & freqs[1] == 0){
                           tmp[, paste("cat", 1:10, sep= "_")] <- NA
                           tmp[, paste("freq", 1:10, sep= "_")] <- NA
                         } else {
                           tmp[, paste("cat", 1:freqs, sep= "_")] <- names(n.categories[[x]])[1:freqs]
                           tmp[, paste("freq", 1:freqs, sep= "_")] <- unclass(n.categories[[x]])[1:freqs] 
                         }
                         tmp
                       })
    
    n.output <- data.frame(do.call("rbind", n.output))
    
    write.csv(n.output, out.file.cat, row.names= TRUE)
    cat('Check for categorical variables completed')
    cat(' // ')
    cat('Results saved to disk') 
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat(' // ')
  print(time.taken)
}