#' Convenience Wrapper to SQL classics for BEGIN,COMMIT, ROLLBACK
#' 
#' this set of function can speed up loops by starting a transaction, 
#' performing several queries and ending them with either commit or rollback. 
#' 
#' @param con PostgreSQL connection object.
#' @rdname transactionUtils
#' @export
beginTransaction <- function(con){
  if(is.null(DBI::dbGetQuery(con,'BEGIN'))) print('BEGIN')
}

#' @rdname transactionUtils
#' @export
commitTransaction <- function(con){
  if(is.null(DBI::dbGetQuery(con,'COMMIT'))) print('COMMIT')
}

#' @rdname transactionUtils
#' @export
rollbackTransaction <- function(con){
  if(is.null(DBI::dbGetQuery(con,'COMMIT'))) print('ROLLBACK')
}
