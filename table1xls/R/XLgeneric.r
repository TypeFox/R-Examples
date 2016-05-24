##' Export a generic data frame, matrix or table to a spreadsheet and save the file.
##'

##' This function is a convenience wrapper for getting practically any rectangular data structure into a spreadsheet, without worrying about conversion or spreadsheet-writing technicalities.
##' 
##' If the structure is not a data frame (or inherited from one), but a table or matrix, the function will convert it into one using \code{\link{as.data.frame.matrix}}, because data frames are what the underlying function \code{\link{writeWorksheet}} can export.
##' 
##' See the \code{\link{XLtwoWay}} help page, for behavior regarding new-sheet creation, overwriting, etc.
##' 
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
##' @seealso For two-way contingency tables, see \code{\link{XLtwoWay}}.
##' 
##' 
##' @example inst/examples/ExGeneric.r
##' 
##' @title Write generic rectangular data to a spreadsheet
##' 
##' @param wb a \code{\link[XLConnect]{workbook-class}} object
##' @param sheet numeric or character: a worksheet name (character) or position (numeric) within \code{wb}. 
##' @param dataset the rectangular structure to be written. Can be a data frame, table, matrix or similar.
##' @param title character: an optional overall title to the table. Default (\code{NULL}) is no title.
##' @param addRownames logical: should a column of row names be added to the left of the structure? (default \code{FALSE})
##' @param rowTitle character: the title to be placed above the row name column (default "Name")
##' @param rowNames character: vector of row names. Default \code{rownames(dataset)}, but relevant only if \code{addRownames=TRUE}.
##' @param colNames character: vector of column names to replace the original ones. Default \code{NULL}, meaning that the original names are left intact. Note that the title for the row-names column (if \code{addRownames=TRUE}) is *not* considered part of \code{colNames}, and is set separately.
##' @param row1,col1 numeric: the first row and column occupied by the output. 
##' @param purge logical: should \code{sheet} be created anew, by first removing the previous copy if it exists? (default \code{FALSE})
##' 
##' @return The function returns invisibly, after writing the data into \code{sheet} and saving the file.
##'
##' @export
##' @import XLConnect

XLgeneric<-function(wb,sheet,dataset,title=NULL,addRownames=FALSE,rowNames=rownames(dataset),rowTitle="Name",colNames=NULL,row1=1,col1=1,purge=FALSE)
{ 
  
  if (!("data.frame" %in% class(dataset))) dataset=as.data.frame.matrix(dataset)

### Optionally changing the column names  
  if(!is.null(colNames))
  {
    if(length(colNames)!=dim(dataset)[2]) stop("Number of column names doesn't match dataset width.\n")
    names(dataset)=colNames
  }  
  
  if(purge) removeSheet(wb,sheet)
  if(!existsSheet(wb,sheet)) createSheet(wb,sheet)
  
  if(!is.null(title))  ### Adding a title
  {
    XLaddText(wb,sheet,text=title,row1=row1,col1=col1)
    row1=row1+1
  }
  if(addRownames)
  {
    if(length(rowNames)!=dim(dataset)[1]) stop("Number of row names doesn't match dataset length.\n")
    dataset=cbind(rowNames,dataset)
      names(dataset)[1]=rowTitle
  }
  writeWorksheet(wb,dataset,sheet,startRow=row1,startCol=col1)
  

setColumnWidth(wb, sheet = sheet, column = col1:(dim(dataset)[2]+1), width=-1)

saveWorkbook(wb)
}  ### Function end
