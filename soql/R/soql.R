soql <- function(query = '') {
  soql_list <- list(clauses = list(), simple_filters = list(), endpoint = NULL)
  class(soql_list) <- 'soql'
  
  endpoint_and_params <- strsplit(query, '?', fixed = TRUE)[[1]]
  if(!is.na(endpoint_and_params[1])) {
    soql_list$endpoint <- endpoint_and_params[1]
  }
  params <- endpoint_and_params[2]
  
  if(is.na(params)) {
    # No need to do the rest
    return(soql_list)
  }
  
  params <- strsplit(params, '&', fixed = TRUE)[[1]]
  
  for(param in params) {
    param_key_value <- strsplit(param, '=', fixed = TRUE)[[1]]
    param_key <- URLdecode(param_key_value[1])
    param_value <- URLdecode(param_key_value[2])
    
    as_numeric <- suppressWarnings({as.numeric(param_value)})
    if(!is.na(as_numeric)) {
      param_value <- as_numeric
    }
    
    if(substring(param_key, 1, 1) == '$') {
      soql_list$clauses[[substring(param_key, 2)]] <- param_value
    } else {
      soql_list$simple_filters[[param_key]] <- param_value
    }
  }
  return(soql_list)
}

as.character.soql <- function(x, ...) {
  params <- '?'
  
  for(param_key in names(x$clauses)) {
    if(params != '?') {
      params <- paste0(params, '&')
    }
    
    param_value <- as.character(x$clauses[[param_key]])
    
    params <- paste0(
                params,
                '$',
                URLencode(param_key),
                '=',
                URLencode(param_value)
              )
  }
  
  for(param_key in names(x$simple_filters)) {
    if(params != '?') {
      params <- paste0(params, '&')
    }
    
    param_value <- as.character(x$simple_filters[[param_key]])
    
    params <- paste0(
      params,
      URLencode(param_key),
      '=',
      URLencode(param_value)
    )
  }
  
  if(params == '?') {
    return(x$endpoint)
  } else {
    # Null resolves to ''
    return(paste0(x$endpoint, params))
  }
}

print.soql <- function(x, ...) {
  result <- as.character(x)
  if(!is.null(result)) {
    cat(result)
  }
}

.concat_clause <- function(clause, new_clause_segment, sep = ',') {
  if(is.null(clause)) {
    return(new_clause_segment)
  } else {
    return(paste0(clause, sep, new_clause_segment))
  }
}

soql_where <- function(soql_list, where_clause) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  # If there's an OR in it, needs to be parenthesized
  if(!is.null(where_clause) && grepl('or', tolower(where_clause))) {
    where_clause <- paste0('(', where_clause, ')')
  }
  
  soql_list$clauses$where <- .concat_clause(soql_list$clauses$where, where_clause, sep = ' AND ')
  
  return(soql_list)
}

soql_select <- function(soql_list, select_clause) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$clauses$select <- .concat_clause(soql_list$clauses$select, select_clause)
  
  return(soql_list)
}

soql_order <- function(soql_list, column, desc = FALSE) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  if(desc) {
    column <- paste(column, 'DESC')
  }

  soql_list$clauses$order <- .concat_clause(soql_list$clauses$order, column)
  
  return(soql_list)
}

soql_group <- function(soql_list, group_clause) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$clauses$group <- .concat_clause(soql_list$clauses$group, group_clause)
  
  return(soql_list)
}

soql_limit <- function(soql_list, limit) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }

  soql_list$clauses$limit <- as.numeric(limit)
  
  return(soql_list)
}

soql_offset <- function(soql_list, offset) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$clauses$offset <- as.numeric(offset)
  
  return(soql_list)
}

soql_q <- function(soql_list, q_clause) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$clauses$q <- .concat_clause(soql_list$clauses$q, q_clause, sep = ' ')
  
  return(soql_list)
}

soql_simple_filter <- function(soql_list, column, value) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$simple_filters[[column]] <- as.character(value)
  
  return(soql_list)
}

soql_add_endpoint <- function(soql_list, endpoint) {
  if(class(soql_list) != 'soql') {
    return(NULL)
  }
  
  soql_list$endpoint <- endpoint
  
  return(soql_list)
}