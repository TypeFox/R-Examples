##' One-way Contingency Tables exported to a spreadsheet
##'
##' Calculates a one-way contingency table in counts and percents, exports a formatted output to a spreadsheet, and saves the file.
##' 
##' This function performs a one-way contingency table, also calculating the distribution in percents. 
##' 
##' The table is then exported to worksheet \code{sheet} in workbook \code{wb}, either using the format \code{"Count(percent)"} (if \code{combine=TRUE}), or as two separate columns in the same table. 
##' 
##' See the \code{\link{XLtwoWay}} help page, for behavior regarding new-sheet creation, overwriting, etc.
##' 
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
##' @seealso If interested in other descriptive summaries, see \code{\link{XLunivariate}}. For two-way contingency tables, see \code{\link{XLtwoWay}}.
##' 
##' 
##' @example inst/examples/Ex1way.r 
##' 
##' 
##' @param wb a \code{\link[XLConnect]{workbook-class}} object
##' @param sheet numeric or character: a worksheet name (character) or position (numeric) within \code{wb}. 
##' @param rowvar vector: the categorical variable (logical, numeric, character, factor, etc.) to be tabulated
##' @param title character: an optional overall title to the table. Default (\code{NULL}) is no title.
##' @param rowTitle character: the title to be placed above the row name column (default empty string)
##' @param rowNames character: vector of row names. Default behavior (\code{NULL}): automatically determined from data
##' @param ord numeric vector specifying row-index order in the produced table. Default (\code{NULL}) is no re-ordering.
##' @param row1,col1 numeric: the first row and column occupied by the table (title included if relevant).
##' @param purge logical should \code{sheet} be created anew, by first removing the previous copy if it exists? (default \code{FALSE})
##' @param digits numeric: how many digits (after the decimal point) to show in the percents? Defaults to 1 if n>=200, 0 otherwise.
##' @param combine logical: should counts and percentab be combined to the popular \code{"Count(percent)"} format, or presented side-by-side? (default \code{TRUE}) 
##' @param useNA How to handle missing values. Passed on to \code{\link{table}} (see help on that function for options).
##' 
##' @return The function returns invisibly, after writing the data into \code{sheet} and saving the file.
##'
##' @export

XLoneWay<-function(wb,sheet,rowvar,title=NULL,rowTitle="Value",rowNames=NULL,ord=NULL,row1=1,col1=1,purge=FALSE,digits=ifelse(length(rowvar)>=200,1,0),combine=TRUE,useNA='ifany')
{ 
  
  if(purge) removeSheet(wb,sheet)
  if(!existsSheet(wb,sheet)) createSheet(wb,sheet)
  
  if(!is.null(title))  ### Adding a title
  {
    XLaddText(wb,sheet,text=title,row1=row1,col1=col1)
    row1=row1+1
  }
    
  n=length(rowvar)
  tab=table(rowvar,useNA=useNA)
  percentab=round(tab*100/n,digits=digits)
  names(tab)[is.na(names(tab))]="missing"
  tabnames=names(tab)
  tab=c(tab,n)
  names(tab)=c(tabnames,"Total")
  percentab=c(percentab,100)
  if (is.null(ord)) ord=1:length(tab)
  if (!is.null(rowNames)) names(tab)=rowNames

if(combine) {
    
    dout=data.frame(cbind(names(tab),paste(tab,' (',percentab,'%)',sep="")))
    names(dout)=c(rowTitle,"Count (%)")
    
} else {

  dout=data.frame(nam=names(tab),Count=tab,Percent=percentab)
  names(dout)[1]=rowTitle
}
  writeWorksheet(wb,dout[ord,],sheet,startRow=row1,startCol=col1)   

setColumnWidth(wb, sheet = sheet, column = col1:(col1+5), width=-1)
saveWorkbook(wb)
}  ### Function end
