
turn_wdgt_on <- function(widget_list)
{
  
  lapply(widget_list , FUN = function(h) {enabled(h) <- T})
  
}

turn_wdgt_off <- function(widget_list)
{
  
  lapply(widget_list , FUN = function(h) {enabled(h) <- F})
  
}

get_wdgt_vals <- function(widget_list)
{
  
  as.character(sapply(widget_list, FUN = svalue))
  
}

set_wdgt_sel <- function(widget_list, col_vals)
{
  
  for(i in 1:length(col_vals))
  {
    svalue(widget_list[[i]]) <- col_vals[i]
  }
  
}

set_wdgt_vals <- function(widget_list, col_vals)
{
  
  sapply(widget_list, FUN = function(h) {h[] <- col_vals})
  
}