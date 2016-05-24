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

#' @title Overlay over \link[RODBC]{odbcFetchRows} 
#' @useDynLib RODBCext
#' @import RODBC
#' @description
#' RODBC::odbcFetchRows crashes if the ODBC channel is in "query prepared
#' but already not executed" state.
#' This function is a small overlay emmitting an error in such a case.
#' @param channel ODBC connection obtained by \link{odbcConnect}
#' @param ... other parametrs passed to \link[RODBC]{odbcFetchRows}
#' @return see \link[RODBC]{odbcFetchRows}
#' @export
odbcFetchRows <- function(channel, ...)
{
  if(!odbcValidChannel(channel)){
    stop("first argument is not an open RODBC channel")
  }
  res <- .Call("RODBCQueryStatus", attr(channel, "handle_ptr"))
  
  if(res == -1){
    stop('an error occured')
  }else if(res == 0){
    stop('query is prepared but not executed, you must execute it first')
  }

  return(RODBC::odbcFetchRows(channel, ...))
}