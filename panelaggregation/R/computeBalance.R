#' Compute Balances from a Item Shares
#' 
#' This function computes balances (i.e. positive - negative items), from item shares
#' stored in a wide format data.table.
#' 
#' @author Matthias Bannert, Gabriel Bucu
#' @param data_table a data.table in wide format containing item 
#' @param multipliers list containing multipliers of items, assigned by item and column names
#' @export
 computeBalance <- function(data_table, multipliers = list("item_pos" = 1, "item_eq" = 0, "item_neg" = -1)) {
  data_table$balance <- 0
  
  # get rid of the CRAN check NOTE, this only for getting the package CRAN ready
  # see Matthew Dowle on
  #http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
  balance = NULL 
  
  for(item in names(multipliers)) {
    data_table[, balance := balance + multipliers[[item]] * get(item)]
  }
  
  data_table
}
