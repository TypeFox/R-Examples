tabmulti <- function(dataset, xvarname, yvarnames, ymeasures = NULL, listwise.deletion = TRUE,
                     latex = FALSE, xlevels = NULL, ynames = yvarnames, ylevels = NULL,
                     quantiles = NULL, quantile.vals = FALSE, parenth.sep = "-", decimals = NULL,
                     cell = "n", freq.parenth = NULL, freq.text.label = NULL, freq.tests = "chi", 
                     means.parenth = "sd", means.text.label = NULL, variance = "unequal", 
                     medians.parenth = "iqr", medians.text.label = NULL, p.include = TRUE, 
                     p.decimals = c(2, 3), p.cuts = 0.01, p.lowerbound = 0.001, p.leading0 = TRUE, 
                     p.avoid1 = FALSE, overall.column = TRUE, n.column = FALSE, n.headings = TRUE, 
                     compress = FALSE, bold.colnames = TRUE, bold.varnames = FALSE, 
                     bold.varlevels = FALSE, variable.colname = "Variable") {
  
  # If any inputs are not correct class, return error
  if (!is.matrix(dataset) & !is.data.frame(dataset)) {
    stop("For dataset input, please enter matrix or data frame with variables of interest")
  }
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
  if (!is.null(quantiles) && ! (is.numeric(quantiles) & length(quantiles) == 1 && quantiles > 1 & quantiles == round(quantiles))) {
    stop("For quantiles input, please enter a whole number greater than 1")
  }
  if (!is.logical(quantile.vals)) {
    stop("For quantile.vals input, please enter TRUE or FALSE")
  }
  if (!is.character(parenth.sep)) {
    stop("For parenth.sep input, please enter a character string")
  }
  if (!is.null(decimals) && !all(is.numeric(decimals))) {
    stop("For decimals input, please enter numeric value or vector of numeric values indicating how many decimal places should
         be used in reporting statistics for each row variable")
  }
  if (! cell %in% c("n", "percent", "tot.percent", "col.percent", "row.percent", "tot.prop",
                    "col.prop", "row.prop", "n/totn", "n/coln", "n/rown")) {
    stop("For cell input, please enter 'n', 'tot.percent', 'col.percent', 'row.percent', 
         'tot.prop', 'col.prop', 'row.prop', 'n/totn', 'n/coln', or 'n/rown'")
  }
  if (!is.null(freq.parenth) && ! freq.parenth %in% c("none", "se", "ci", "tot.percent", "col.percent", "row.percent", 
                                                      "tot.prop", "col.prop", "row.prop")) {
    stop("For freq.parenth input, please enter 'none', 'se', 'ci', 'tot.percent', 'col.percent', 'row.percent', 
         'tot.prop', 'col.prop', or 'row.prop'")
  }
  if (!is.null(freq.text.label) && !is.character(freq.text.label)) {
    stop("For freq.text.label input, please enter a character string or just leave it unspecified. Use 'none' to request no label")
  }
  if (!all(freq.tests %in% c("chi", "fisher", "z", "z.continuity"))) {
    stop("For freq.tests input, please enter character string or vector of character strings indicating what statistical test
         should be performed for each categorical row variable. Each element should be 'chi', 'fisher', 'z', or 'z.continuity'")
  }
  if (! means.parenth %in% c("none", "sd", "se", "t.ci", "z.ci", "none")) {
    stop("For means.parenth input, please enter 'none', 'sd', 'se', 't.ci', or 'z.ci'")
  }
  if (!is.null(means.text.label) && !is.character(means.text.label)) {
    stop("For means.text.label input, please enter a character string or just leave it unspecified. Use 'none' to request no label")
  }
  if (! variance %in% c("equal", "unequal", "ftest")) {
    stop("For variance input, please enter 'equal', 'unequal', or 'ftest'")
  }
  if (! medians.parenth %in% c("none", "iqr", "range", "minmax", "q1q3")) {
    stop("For medians.parenth input, please enter 'none', 'iqr', 'range', 'minmax', or 'q1q3'")
  }
  if (!is.null(medians.text.label) && !is.character(medians.text.label)) {
    stop("For medians.text.label input, please enter a character string or just leave it unspecified. Use 'none' to request no label")
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
  if (!is.logical(overall.column)) {
    stop("For overall.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.column)) {
    stop("For n.column input, please enter TRUE or FALSE")
  }
  if (!is.logical(n.headings)) {
    stop("For n.headings input, please enter TRUE or FALSE")
  }
  if (!is.logical(compress)) {
    stop("For compress input, please enter TRUE or FALSE")
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
  
  # If listwise.deletion is TRUE, drop observations with missing values for column variable or any row variables
  if (listwise.deletion == TRUE){
    
    dataset <- dataset[which(!is.na(dataset[, xvarname])), ]
    for (ii in 1:length(yvarnames)) {
      dataset <- dataset[which(!is.na(dataset[, yvarnames[ii]])), ]
    }
    
  }
  
  # If ymeasures is NULL, guess what measures are appropriate based on each variable
  if (is.null(ymeasures)) {
    ymeasures <- c()
    for (ii in 1:length(yvarnames)) {
      if (is.factor(dataset[, yvarnames[ii]]) | length(unique(dataset[!is.na(dataset[, yvarnames[ii]]), yvarnames[ii]])) <= 5) {
        ymeasures <- c(ymeasures, "freq")
      } else {
        ymeasures <- c(ymeasures, "mean")
      }
    }
  }
  
  # If ymeasures is single value, create vector of repeat values
  if (length(ymeasures) == 1) {
    ymeasures <- rep(ymeasures, length(yvarnames))
  }
  
  # If freq.tests is a single value, create vector of repeat values
  if (length(freq.tests) == 1) {
    freq.tests <- rep(freq.tests, sum(ymeasures == "freq"))
  }
  
  # If decimals is a single value, create vector of repeat values
  if (length(decimals) == 1) {
    decimals <- rep(decimals, length(ymeasures))
  }
  
  # If ylevels is a vector, convert to a list
  if (!is.null(ylevels) && !is.list(ylevels)) {
    ylevels <- list(ylevels)
  }
  
  # Call tabmeans, tabmedians, or tabfreq repeatedly
  mediansindex <- 0  
  meansindex <- 0  
  freqindex <- 0
  for (ii in 1:length(yvarnames)) {
    if (ymeasures[ii] == "mean") {
      meansindex <- meansindex + 1
      current <- tabmeans(x = dataset[, xvarname], y = dataset[, yvarnames[ii]], latex = latex, variance = variance, 
                          xlevels = xlevels, yname = ynames[ii], quantiles = quantiles, quantile.vals = quantile.vals, 
                          parenth = means.parenth, text.label = means.text.label, parenth.sep = parenth.sep,
                          decimals = decimals[ii], p.include = p.include, p.decimals = p.decimals, p.cuts = p.cuts,
                          p.lowerbound = p.lowerbound, p.leading0 = p.leading0, p.avoid1 = p.avoid1, 
                          overall.column = overall.column, n.column = n.column, n.headings = n.headings, 
                          bold.colnames = bold.colnames, bold.varnames = bold.varnames, variable.colname = variable.colname)
    } else if (ymeasures[ii] == "median") {
      mediansindex <- mediansindex + 1
      current <- tabmedians(x = dataset[, xvarname], y = dataset[, yvarnames[ii]], latex = latex, xlevels = xlevels,
                            yname = ynames[ii], quantiles = quantiles, quantile.vals = quantile.vals, 
                            parenth = medians.parenth, text.label = medians.text.label, parenth.sep = parenth.sep, 
                            decimals = decimals[ii], p.include = p.include, p.decimals = p.decimals, p.cuts = p.cuts, 
                            p.lowerbound = p.lowerbound, p.leading0 = p.leading0, p.avoid1 = p.avoid1, 
                            overall.column = overall.column, n.column = n.column, n.headings = n.headings, 
                            bold.colnames = bold.colnames, bold.varnames = bold.varnames, variable.colname = variable.colname)
    } else if (ymeasures[ii] == "freq") {
      freqindex <- freqindex + 1
      current <- tabfreq(x = dataset[, xvarname], y = dataset[, yvarnames[ii]], latex = latex, xlevels = xlevels,
                         yname = ynames[ii], ylevels = ylevels[[freqindex]], quantiles = quantiles, 
                         quantile.vals = quantile.vals, cell = cell, parenth = freq.parenth, text.label = freq.text.label,
                         parenth.sep = parenth.sep, test = freq.tests[freqindex], decimals = decimals[ii], p.include = p.include,
                         p.decimals = p.decimals, p.cuts = p.cuts, p.lowerbound = p.lowerbound, p.leading0 = p.leading0, 
                         p.avoid1 = p.avoid1, overall.column = overall.column, n.column = n.column, n.headings = n.headings, 
                         compress = compress, bold.colnames = bold.colnames, bold.varnames = bold.varnames, 
                         bold.varlevels = bold.varlevels, variable.colname = variable.colname)
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