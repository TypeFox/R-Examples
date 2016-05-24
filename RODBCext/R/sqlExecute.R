# Copyright (C) 2014 Mateusz Zoltak
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

#' @title Executes an already prepared query
#' @useDynLib RODBCext
#' @import RODBC
#' @description
#' Executes a parameterized query. 
#' 
#' Optionally (fetch=TRUE) fetches results using \link[RODBC]{sqlGetResults}.
#' 
#' Optionally (query=NULL) uses query already prepared by \link{sqlPrepare}.
#' @details
#' Return value depends on the combination of parameters:
#' \itemize{
#'   \item if there were errors during query preparation or execution or fetching results
#'     return value depends on errors parameter - if errors=TRUE error is thrown,
#'     otherwise -1 will be returned
#'   \item if fetch=FALSE and there were no errors invisible(1) will be returned
#'   \item if fetch=TRUE and there were no errors a data.frame with results will be returned
#' }
#' @param channel ODBC connection obtained by \link[RODBC]{odbcConnect}
#' @param query a query string (NULL if query already prepared using
#'   \link{sqlPrepare})
#' @param data data to pass to sqlExecute (as data.frame)
#' @param fetch whether to automatically fetch results (if data provided)
#' @param errors whether to display errors
#' @param rows_at_time number of rows to fetch at one time - see details of
#'   \link[RODBC]{sqlQuery}
#' @param force_loop whether to execute queries in the explicit loop with 
#'   separate query planing for each iteration (usefull if executing a query 
#'   invalidates its plan, e.g. EXEC queries on Ms SQL Server)
#' @param ... parameters to pass to \link[RODBC]{sqlGetResults} (if fetch=TRUE)
#' @return see datails
#' @export
#' @examples
#' \dontrun{
#'   conn = odbcConnect('MyDataSource')
#'   
#'   # prepare, execute and fetch results separatly
#'   sqlPrepare(conn, "SELECT * FROM myTable WHERE column = ?")
#'   sqlExecute(conn, NULL, 'myValue')
#'   sqlGetResults(conn)
#'   
#'   # prepare and execute at one time, fetch results separately
#'   sqlExecute(conn, "SELECT * FROM myTable WHERE column = ?", 'myValue')
#'   sqlGetResults(conn)
#'   
#'   # prepare, execute and fetch at one time
#'   sqlExecute(conn, "SELECT * FROM myTable WHERE column = ?", 'myValue', TRUE)
#'   
#'   # prepare, execute and fetch at one time, pass additional parameters to sqlFetch()
#'   sqlExecute(
#'     conn, 
#'     "SELECT * FROM myTable WHERE column = ?", 
#'     'myValue', 
#'     TRUE, 
#'     stringsAsFactors=FALSE
#'   )
#' }
sqlExecute <- function(channel, query=NULL, data=NULL, fetch=FALSE, errors=TRUE, rows_at_time=attr(channel, "rows_at_time"), force_loop=FALSE, ...)
{
  stopifnot(
    odbcValidChannel(channel),
    is.vector(fetch), is.logical(fetch), length(fetch) == 1, all(!is.na(fetch)),
    is.vector(errors), is.logical(errors), length(errors) == 1, all(!is.na(errors)),
    is.vector(rows_at_time), is.numeric(rows_at_time), length(rows_at_time) == 1, all(!is.na(rows_at_time)),
    is.vector(force_loop), is.logical(force_loop), length(force_loop) == 1, all(!is.na(force_loop))
  )
  if(!odbcValidChannel(channel)){
    stop("first argument is not an open RODBC channel")
  }
  
  # workaround for queries wchich has to be planned before each execution
  if(force_loop){
    data = as.data.frame(data)
    stopifnot(
      is.vector(query), is.character(query), length(query) == 1, all(!is.na(query)),
      nrow(data) > 0
    )
    results = list()
    for(i in seq_along(data[, 1])){
      results[[i]] = sqlExecute(channel, query, data[i, ], fetch, errors, rows_at_time, FALSE, ...)
    }
    return(do.call(rbind, results))
  }
  
  # Prepare query (if provided)
  if(!is.null(query)){
    stat <- sqlPrepare(channel, query, errors)
    if(stat == -1L){
      return(stat); # there is no need to check if error should be thrown - this is being done by sqlPrepare()
    }
  }
  
  # Prepare data
  data = as.data.frame(data)
  for(k in seq_along(data)){
    if(is.factor(data[, k])){
      data[, k] = levels(data[, k])[data[, k]]
    }
    if(is.logical(data[, k])){
      data[, k] = as.character(data[, k])
    }
  }
  
  # If there is no need to fetch results or no query parameters were provided,
  # call RODBCExecute once on whole data
  if(fetch == FALSE | nrow(data) < 1){
    stat <- .Call(
      "RODBCExecute", 
      attr(channel, "handle_ptr"), 
      data, 
      as.integer(rows_at_time)
    )
    if(stat == -1L) {
      if(errors){
        stop(paste0(RODBC::odbcGetErrMsg(channel), collapse='\n'))
      }
      else{
        return(stat)
      }
    }
    
    if(fetch == FALSE){
      return(invisible(stat))
    }
    
    # Fetch results
    return(RODBC::sqlGetResults(channel, errors=errors, ...))
  }
  
  # If results should be fetched and query parameters were provided

  # For each row of query parameters execute the query and fetch results
  results = NULL
  for(row in 1:nrow(data)){
    stat <- .Call(
      "RODBCExecute", 
      attr(channel, "handle_ptr"), 
      as.list(data[row, ]), 
      as.integer(rows_at_time)
    )
    if(stat == -1L) {
      if(errors){
        stop(paste0(RODBC::odbcGetErrMsg(channel), collapse='\n'))
      }
      else{
        return(stat)
      }
    }      
    
    stat <- RODBC::sqlGetResults(channel, errors=errors, ...)
    
    if(is.null(results)){
      results <- as.data.frame(stat)
    }else{
      results <- rbind(results, as.data.frame(stat))
    }
  }
  return(results)
}
