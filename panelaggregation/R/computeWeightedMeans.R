#' Compute Weighted Mean by Group
#' 
#' This function computes the weighted mean of variable groups from a data.table. 
#' computeWeightedMean is performance optimized and designed to
#' work well in bulk operations. The function returns a data.table.
#' 
#' @author Matthias Bannert, Gabriel Bucur
#' @param data_table a data.table
#' @param variables character name of the variable(s) to focus on. The variables must be in the data.table
#' @param weight character name of the data.table column that contains a weight. 
#' @param by character vector of the columns to group by
#' @import data.table
#' @example demo/aggregation.R
#' @export
computeWeightedMeans <- function(data_table, variables, weight, by) {
  
  if (is.null(weight)) {
    res_dt <- data_table[, lapply(.SD, mean, na.rm = TRUE), .SDcols = variables, by = by]
  } else {
    res_dt <- data_table[, lapply(.SD, weighted.mean, weight = eval(as.name(weight)), na.rm = TRUE),
                         .SDcols = variables, by = by]
  }
  
  res_dt
}
