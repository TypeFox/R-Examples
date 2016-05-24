tabmulti.svy <- function(svy, xvarname, yvarnames, ymeasures = NULL, listwise.deletion = FALSE,
                         latex = FALSE, xlevels = NULL, ynames = yvarnames, ylevels = NULL, 
                         mean.tests = "Wald", median.tests = "wilcoxon", freq.tests = "F", 
                         decimals = 1, p.include = TRUE, p.decimals = c(2, 3), p.cuts = 0.01, 
                         p.lowerbound = 0.001, p.leading0 = TRUE, p.avoid1 = FALSE, n.column = FALSE, 
                         n.headings = TRUE, se = FALSE, compress = FALSE, parenth = "iqr", 
                         text.label = NULL, parenth.sep = "-", bold.colnames = TRUE, 
                         bold.varnames = FALSE, bold.varlevels = FALSE, variable.colname = "Variable") {
  
  # If any inputs are not correct class, return error
  if (!is.character(xvarname)) {
    stop("For xvarname input, please enter character string with name of column variable")
  }
  if (!all(is.character(yvarnames))) {
    stop("For yvarnames input, please enter character string or vector of character strings with name(s) of row variable(s)")
  }
  if (!is.null(ymeasures) && !all(ymeasures %in% c("mean", "median", "freq"))) {
    stop("For ymeasures input, please enter character string or vector of character strings of same length as yvarnames")
  }
  if (!is.logical(listwise.deletion)) {
    stop("For listwise.deletion input, please enter TRUE or FALSE")
  }
  if (!is.logical(latex)) {
    stop("For latex input, please enter TRUE or FALSE")
  }
  if (!is.null(xlevels) && !is.character(xlevels)) {
    stop("For xlevels input, please enter vector of character strings")
  }
  if (!all(is.character(ynames))) {
    stop("For ynames input, please enter character string or vector of character strings of same length as yvarnames")
  }
  if (!is.null(ylevels) && !all(unlist(lapply(X = ylevels, FUN = function(x) all(is.character(x)))))) {
    stop("For ylevels input, please enter vector or list of vectors of character strings")
  }
  if (!all(mean.tests %in% c("Wald", "LRT"))) {
    stop("For mean.tests input, please enter character string or vector of character strings specifying whether Wald or 
          Likelihood Ratio Test statistic should be used for each mean comparison. Each element should be 'Wald' or 'LRT'.")
  }
  if (!all(median.tests %in% c("wilcoxon", "vanderWaerden", "median", "KruskalWallis"))) {
    stop("For median.tests input, please enter character string or vector of character strings indicating what statistical
         test should be used for each median comparison. Each element should be a possible value for the 'test' input of the 
         svyranktest function in the survey package: 'wilcoxon', 'vanderWaerden', 'median', or 'KruskalWallis'. See 
         documentation for tabmedians.svy and svyranktest for details.")
  }     
  if (!all(freq.tests %in% c("F", "Chisq", "Wald", "adjWald", "lincom", "saddlepoint"))) {
    stop("For freq.tests input, please enter character string or vector of character strings indicating what statistical test
         should be performed for each categorical row variable. Each element should be a possible value for the 'statistic' 
         input of the svychisq function in the survey package: 'F', 'Chisq', 'Wald', 'adjWald', 'lincom', 
         or 'saddlepoint'. See svychisq documentation for details.")
  }
  if (!is.numeric(decimals)) {
    stop("For decimals input, please enter numeric value")
  }
  if (!is.logical(p.include)) {
    stop("For p.include input, please enter TRUE or FALSE")
  }
  if (!is.numeric(p.decimals)) {
    stop("For p.decimals input, please enter numeric value or vector")
  }
  if (!is.numeric(p.cuts)) {  
    stop("For p.cuts input, please enter numeric value or vector")
  }
  if (!is.numeric(p.lowerbound)) {
    stop("For p.lowerbound input, please enter numeric value")
  }
  if (!is.logical(p.leading0)) {
    stop("For p.leading0 input, please enter TRUE or FALSE")
  }
  if (!is.logical(p.avoid1)) {
    stop("For p.avoid1 input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.column)) {
    stop("For n.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.headings)) {
    stop("For n.headings input, please enter TRUE or FALSE")
  }
  if (!is.logical(se)) {
    stop("For se input, please enter TRUE or FALSE")
  }
  if (!is.logical(compress)) {
    stop("For compress input, please enter TRUE or FALSE")
  }
  if (! parenth %in% c("minmax", "range", "q1q3", "iqr", "none")) {
    stop("For parenth input, please enter one of the following: 'minmax', 'range', 'q1q3', 'iqr', 'none'")
  }
  if (!is.null(text.label) && !is.character(text.label)) {
    stop("For text.label input, please enter something like 'Median (IQR)' or just leave it unspecified")
  }
  if (!is.character(parenth.sep)) {
    stop("For parenth.sep input, please enter a character string")
  }
  if (!is.logical(bold.colnames)) {
    stop("For bold.colnames input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.varnames)) {
    stop("For bold.varnames input, please enter TRUE or FALSE")
  }
  if (!is.logical(bold.varlevels)) {
    stop("For bold.varlevels input, please enter TRUE or FALSE")
  }
  if (!is.character(variable.colname)) {
    stop("For variable.colname input, please enter a character string")
  }
  
  # Save xvarname and yvarnames character strings
  xstring <- xvarname
  x <- svy$variables[, xstring]
  
  # If listwise.deletion is TRUE, drop observations with missing values for column variable or any row variables
  if (listwise.deletion == TRUE) {
    
    # Loop through and find locations of complete data
    locs <- rep(1, nrow(svy))
    locs[is.na(x)] <- 0
    for (ii in 1:length(yvarnames)) {
      ystring <- yvarnames[ii]
      y <- svy$variables[, ystring]
      locs[is.na(y)] <- 0
    }
    svy <- subset(svy, locs)
    
  }
  
  # If ymeasures is single value, create vector of repeat values
  if (length(ymeasures) == 1) {
    ymeasures <- rep(ymeasures, length(yvarnames))
  }
  
  # If freq.tests is a single value, create vector of repeat values
  if (length(freq.tests) == 1) {
    freq.tests <- rep(freq.tests, length(yvarnames))
  }
  
  # If mean.tests is a single value, create vector of repeat values
  if (length(mean.tests) == 1) {
    mean.tests <- rep(mean.tests, length(yvarnames))
  }
  
  # If median.tests is a single value, create vector of repeat values
  if (length(median.tests) == 1) {
    median.tests <- rep(median.tests, length(yvarnames))
  }
  
  # If ymeasures is NULL, guess what measures are appropriate based on each variable
  if (is.null(ymeasures)) {
    ymeasures <- c()
    for (ii in 1:length(yvarnames)) {
      
      # Save x and y as character strings
      xstring <- xvarname
      ystring <- yvarnames[ii]
      
      # Extract vectors x and y
      x <- svy$variables[, xstring]
      y <- svy$variables[, ystring]
      
      # Find indices for non-missing x and y
      locs <- which(!is.na(x) & !is.na(y))
      x <- x[locs]
      y <- y[locs]
      
      if (is.factor(y) | length(unique(y)) <= 5) {
        ymeasures <- c(ymeasures, "freq")
      } else {
        ymeasures <- c(ymeasures, "mean")
      }
    }
  }
  
  # If ylevels is a vector, convert to a list
  if (!is.null(ylevels) && !is.list(ylevels)) {
    ylevels <- list(ylevels)
  }
  
  # Call tabmeans.svy, tabmedians.svy, or tabfreq.svy repeatedly
  freqindex <- 0
  meanindex <- 0
  medianindex <- 0
  for (ii in 1:length(yvarnames)) {
    if (ymeasures[ii] == "mean") {
      meanindex <- meanindex + 1
      current <- tabmeans.svy(x = xvarname, y = yvarnames[ii], svy = svy, latex = latex, xlevels = xlevels, 
                              yname = ynames[ii], test = mean.tests[meanindex], decimals = decimals, p.decimals = p.decimals, 
                              p.cuts = p.cuts, p.lowerbound = p.lowerbound, p.leading0 = p.leading0, p.avoid1 = p.avoid1, 
                              n.column = n.column, n.headings = n.headings, bold.colnames = bold.colnames, 
                              bold.varnames = bold.varnames, variable.colname = variable.colname)
    } else if (ymeasures[ii] == "median") {
      medianindex <- medianindex + 1
      current <- tabmedians.svy(x = xvarname, y = yvarnames[ii], svy = svy, latex = latex, xlevels = xlevels, 
                                yname = ynames[ii], test = median.tests[medianindex], decimals = decimals, p.decimals = p.decimals, 
                                p.cuts = p.cuts, p.lowerbound = p.lowerbound, p.leading0 = p.leading0, p.avoid1 = p.avoid1, 
                                n.column = n.column, n.headings = n.headings, parenth = parenth, text.label = text.label,
                                parenth.sep = parenth.sep, bold.colnames = bold.colnames, bold.varnames = bold.varnames, 
                                variable.colname = variable.colname)
    } else if (ymeasures[ii] == "freq") {
      freqindex <- freqindex + 1
      current <- tabfreq.svy(x = xvarname, y = yvarnames[ii], svy = svy, latex = latex, xlevels = xlevels, 
                             yname = ynames[ii], ylevels = ylevels[[freqindex]], test = freq.tests[freqindex], 
                             decimals = decimals, p.decimals = p.decimals, p.cuts = p.cuts, 
                             p.lowerbound = p.lowerbound, p.leading0 = p.leading0, p.avoid1 = p.avoid1, 
                             n.column = n.column, n.headings = n.headings, compress = compress, 
                             bold.colnames = bold.colnames, bold.varnames = bold.varnames, 
                             variable.colname = variable.colname)
    }
    
    if (ii == 1) {
      results <- current
    } else {
      results <- rbind(results, current)
    }
  }
  rownames(results) <- NULL
  
  # Return results matrix
  return(results)
  
}