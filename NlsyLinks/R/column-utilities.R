#' @name ColumnUtilities
#' @aliases VerifyColumnExists RenameColumn RenameNlsyColumn
#' @export VerifyColumnExists RenameColumn RenameNlsyColumn

#' @title  A collection of functions that helps data management \code{data.frame}s, particularly those derived from NLSY Extracts.
#' @description  A collection of functions that helps data management \code{data.frame}s, particularly those derived from NLSY Extracts.
#' 

#' @usage VerifyColumnExists( dataFrame, columnName )
#' 
#' RenameColumn( dataFrame, oldColumnName, newColumnName )
#' 
#' RenameNlsyColumn( dataFrame, nlsyRNumber, newColumnName )
#' 
#' @param dataFrame The \code{data.frame} whose columns are to be verified or renamed.
#' @param columnName The name of the column to verify is present in the \code{data.frame}.
#' @param nlsyRNumber The name of the column to change. 
#' @param oldColumnName The name of the column to change. 
#' @param newColumnName The desired name of the column.
#' 
#' @details The RNumber assigned by the NLS has a pattern.  In the Nlsy79 Gen1 dataset, the names start with a `R' or `T' and are followed by seven digits (eg, R0000100).  In the Nlsy79 Gen2 dataset, the names start with `C' or `Y' and are followed by seven digits (eg, C0007030, Y1994600).
#' 
#' In the NLS Investigator, a decimal is present in the RNumber (eg, R00001.00).  When the Investigator saves the dataset as a CSV, the decimal is removed (eg, R0000100).
#' @return *IMPORTANT* The \code{RenameColumn} and \code{RenameNlsyColumn} functions do not use side-effects to rename the \code{data.frame}.  Instead, it returns a new \code{data.frame}.  In the example below, notice the assignment to \code{ds}: \code{ds <- RenameNlsyColumn(...)}.
#' 
#' The \code{VerifyColumnExists} function check that exactly one column exists with the specified \code{columnName}.  If so, the index of the column is returned.  If not, an exception is thrown.
#' 
#' @author Will Beasley
#' 
#' 


VerifyColumnExists <- function( dataFrame, columnName ) {
  indices <- base::match(columnName, colnames(dataFrame))
  if( length(indices) != 1 ) base::stop(paste("Exactly 1 matching column name should be found, but", length(indices), "were found."))
  return( indices )
}

RenameColumn <- function( dataFrame, oldColumnName, newColumnName ) {
  index <- VerifyColumnExists(dataFrame=dataFrame, columnName=oldColumnName)
  base::colnames(dataFrame)[index] <- newColumnName
  return( dataFrame )
}

RenameNlsyColumn <- function( dataFrame, nlsyRNumber, newColumnName ) {
  return( RenameColumn(dataFrame=dataFrame,  oldColumnName=nlsyRNumber, newColumnName=newColumnName) )
}
