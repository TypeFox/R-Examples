# Last modified: 2014-1-25 by Jordan Gressel
#--------------------------------------------

# Output a summary table of all numeric and factor columns,
# with optional output for confidence intervals
DataframeSummary <- function(x, conf.intervals=TRUE) {
  if (class(x) != "data.frame")
    stop("must supply a dataframe")
  
  #-------------------------------------
  # Continuous, numeric data summaries
  x.continuous <- x[ , sapply(x, is.numeric)]
  # If there is only one variable of type numeric, the following code puts it into a dataframe
  if (class(x.continuous)=='numeric') {
    x.continuous <- data.frame(x.continuous)
  }
  if (ncol(x.continuous) > 0) {
    x.min <- apply(x.continuous, 2, min, na.rm=TRUE)
    x.max <- apply(x.continuous, 2, max, na.rm=TRUE)
    x.sd <- apply(x.continuous, 2, sd, na.rm=TRUE)
    x.mean <- as.numeric(lapply(x.continuous, mean, na.rm=TRUE))
    x.n <- apply(!is.na(x.continuous), 2, sum)
    x.na <- apply(is.na(x.continuous), 2, sum)
    row.names <- colnames(x.continuous)
    
    # Build output dataframe
    results.continuous <- data.frame(x.mean, x.sd, x.min, x.max, x.n, x.na)
    colnames(results.continuous) <- c("Mean", "Std Dev", "Min", "Max", "N", "NAs")
    rownames(results.continuous) <- row.names
    
    # Optionally add confidence interval columns 
    if (conf.intervals==TRUE) { 
      calculate.confidence.intervals <- function(x, ...) {  
        error <- qt(0.975,df=length(na.omit(x))-1)*sd(na.omit(x))/sqrt(length(na.omit(x)))
        left <- mean(na.omit(x))-error
        right <- mean(na.omit(x))+error
        c(left, right)
      }
      
      conf.table <- t(apply(x.continuous, 2, calculate.confidence.intervals))
      colnames(conf.table) <- c("95% Conf. (lower)", "95% Conf. (upper)")
      
      results.continuous <- cbind(results.continuous, conf.table)
    }
    
    
    # Print the results
    writeLines("Summary of continuous variables")
    writeLines("-------------------------------")
    print(results.continuous)
  }
  
  #------------------------
  # Categorical summaries
  x.factors <- x[ , sapply(x, is.factor)]
  #If there is only one variable of type factor, the following code puts it into a dataframe
  if (class(x.factors)=='factor') {
    x.factors <- data.frame(x.factors)
  }
  if (ncol(x.factors) > 0) {
    buildCountTable <- function(x) {
      
      Counts <- table(x, useNA="ifany")
      
      # Calculate the percent total for non-missing values
      if (length(x[x==TRUE]) != 0) {
        counts.sans.na <- table(x)
        Percents <- round(100*counts.sans.na/sum(counts.sans.na), 2)
        Percents <- cbind(t(Percents), "<NA>"=0)
      } else {
        Percents <- round(100*Counts/sum(Counts), 2)
      }
      
      output <- rbind(Counts, Percents)
      rownames(output) <- c("Counts", "Percents")
      output
    }
    
    lprint <- function(lst) {
      for (i in 1:length(lst)) {
        cat(names(lst[i]), "\n")
        print(lst[[i]])
        cat("\n")
      }
    }
    
    x.counts <- lapply(x.factors, buildCountTable)
    
    if (ncol(x.continuous) > 0)
      writeLines("\n")
    
    writeLines("Summary of categorical variables")
    writeLines("(first row = categories; second row = counts)")
    writeLines("---------------------------------------------")
    lprint(x.counts)
  }
}