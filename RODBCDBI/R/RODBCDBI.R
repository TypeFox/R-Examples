#' RODBCDBI
#' 
#' Provides Access to Databases Through the ODBC Interface An implementation of R's DBI interface using ODBC package as a back-end. This allows R to connect to any DBMS that has a ODBC driver.
#'
#' @name RODBCDBI
#' @docType package
#' @import methods DBI RODBC
NULL
setOldClass("RODBC")