summaryStats = function(x, ...){
  UseMethod("summaryStats")
}

summaryStats.default = function(x, group = rep("Data", length(x)), 
                               data.order = TRUE, digits = 2, 
                               ...){
  if (!is.factor(group)){
    group = factor(group, 
                   levels = if(data.order){
                              unique(group)
                            }else{
                              sort(unique(group))
                            })
    
  }
  
  getStats = function(x, ...){
    quantiles = quantile(x, probs = c(0.25, 0.5, 0.75), ...)
    names(quantiles) = NULL ## strip off the percentages
    return(c(min = min(x, ...), max = max(x, ...), mean = mean(x, ...),
                var = var(x, ...), sd = sd(x, ...), n = length(x), 
                iqr = IQR(x, ...), skewness = skewness(x, ...),
                lq = quantiles[1], median = quantiles[2], uq = quantiles[3]))
  }
  
  stats = NULL

  if(length(unique(group)) == 1){
    stats = as.list(getStats(x, ...))
    
    with(stats, {
      cat(paste("Minimum value:          " ,round(min, digits), "\n"));
      cat(paste("Maximum value:          " ,round(max, digits), "\n"));
      cat(paste("Mean value:             " ,round(mean, digits), "\n"));
      cat(paste("Median:                 " ,round(median, digits), "\n"));
      cat(paste("Upper quartile:         " ,round(uq, digits), "\n"));
      cat(paste("Lower quartile:         " ,round(lq, digits), "\n"));
      cat(paste("Variance:               " ,round(var, digits), "\n"));
      cat(paste("Standard deviation:     " ,round(sd, digits), "\n"));
      cat(paste("Midspread (IQR):        " ,round(iqr, digits), "\n"));
      cat(paste("Skewness:               " ,round(skewness, digits), "\n"));
      cat(paste("Number of data values:  " ,n, "\n"));
    })
  }else{
    temp.df = data.frame(x = x, group = group)
    stats = aggregate(x~group, data = temp.df, FUN = getStats)
    
    resTable = aggregate(x~group, data = temp.df, FUN = getStats)
    rownames(resTable$x) = resTable$group
    resTable$x = resTable$x[,c("n", "mean", "median", "sd", "iqr")]
    colnames(resTable$x) = c("Sample Size", "Mean", "Median", "Std Dev", "Midspread")
    print(resTable$x)
    
    rownames(stats$x) = stats$group
    stats = as.data.frame(stats$x)
    
  }
  
  invisible(stats)
}

summaryStats.formula = function(x, data = NULL, data.order = TRUE, digits = 2, ...){
  if (missing(x) || (length(x) != 3)){
    stop("missing or incorrect formula")
  }
  
  if(is.null(data)){
    vars = eval(attr(terms(x), "variables"), parent.frame())
  }else{
    vars = eval (attr (terms (x), "variables"), data)
  }
  
  summaryStats(vars[[1]], vars[[2]], data.order = data.order, 
                        digits = digits, ...)
}

summaryStats.matrix = function(x, data.order = TRUE, digits = 2, ...){
  nrows = nrow(x)
  ncols = ncol(x)
  
  x = as.vector(x)
  group = factor(rep(1:ncols, rep(nrows, ncols)))
  
  summaryStats(x, group, data.order = data.order, 
               digits = digits, ...)
}
