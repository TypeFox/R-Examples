#' Compute Weighted Shares By Group
#' 
#' This function computes weighted shares from a data.table. 
#' computeShares is performance optimized and designed to
#' work well in bulk operations. The function returns a data.table.
#' 
#' @author Matthias Bannert, Gabriel Bucur, Oliver Mueller
#' @param data_table a data.table
#' @param variable character name of the variable to focus on. The variable must be in the data.table
#' @param weight character name of the data.table column that contains a weight. 
#' @param by character vector of the columns to group by
#' @param wide logical if true the result is returned in wide format dcast.
#' @import data.table
#' @example demo/aggregation.R
#' @export
computeShares <- function(data_table, variable, weight, by, wide = T) {
  
  old_key <- key(data_table)
  setkeyv(data_table, c(by, variable))
  
  # get rid of the CRAN check NOTE, this only for getting the package CRAN ready
  # see Matthew Dowle on
  #http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  .EACHI = NULL
  
  if (is.null(weight)) {
    # make sure .N is a double because data.table gives it to C which throws an error if integer and double are used
    # in the same division.
    res_dt <- data_table[doUniqueCJ(data_table, c(by, variable)), list(share = as.double(.N)), by = .EACHI]
  } else {
    # do not use get instead of eval(as.name()) here because, if column name equals parameter name (in this case weight)
    # you'll run into a name clash. 
    res_dt <- data_table[doUniqueCJ(data_table, c(by, variable)), list(share = sum(eval(as.name(weight)))), by = .EACHI]
  }
  
  # get rid of the CRAN check NOTE, this only for getting the package CRAN ready
  # see Matthew Dowle on
  #http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  share = NULL
  
  res_dt[is.na(share), share := 0][, share := share / sum(share), by = by]
  

  # we use the wide format by default as functions along the workflow
  # make use of it.
  if(wide){
    f <- as.formula(paste(paste(by, collapse = "+"), "~", variable))
    res_dt <- as.data.table(dcast.data.table(res_dt, f, value.var = "share"))
    possible_answers <- setdiff(names(res_dt), by)
    setnames(res_dt, c(by, paste0("item_", possible_answers)))
  }
  
  setkeyv(data_table, old_key)
  res_dt
}
