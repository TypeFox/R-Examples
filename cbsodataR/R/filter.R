# creates a string that can be used in a $filter query
column_filter <- function(key, values){
  if (is.character(values)){
    values <- paste0("'", values, "'")
  }
  query <- paste0(key, " eq ", values)
  query <- paste0(query, collapse=" or ")
  query
}

# creates a filter string that can be used in $filter query
# usage: get_filter(Periode=c(2001))
get_filter <- function(..., filter_list=list(...)){
  if (length(filter_list) == 0){
    return(NULL)
  }
  
  query <- sapply(names(filter_list), function(key){
    filter <- column_filter(key, filter_list[[key]])
    paste0("(", filter, ")")
  })
  
  paste0(query, collapse=" and ")
}

get_select <- function(select){
  query = NULL
  if (length(select) > 0){
    query = paste0(select, collapse = ",")
  }
  query
}

get_query <- function(..., select=NULL){
  query <- ""
  filter <- get_filter(...)
  if (!is.null(filter)){
    query = paste0(query, "&$filter=", filter)
  }
  
  select <- get_select(select)
  if (!is.null(select)){
    query = paste0(query, "&$select=", select)
  }
  query
}


#column_filter("Periode", c(2000, 2001))
#get_filter(list(Periode=2001:2003), id="2")
#get_query(list(Periode=2001:2003))

