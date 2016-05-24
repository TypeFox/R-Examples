


# This function is internal
fancytab2<-function(x,y=NULL,digits=1,sumby=2,rowvar="",rowNames=NULL,missings='ifany')
{
tout=addmargins(table(x,y,useNA=missings))
pout=round(200*prop.table(tout,margin=sumby),digits)
rownames(tout)[is.na(rownames(tout))]="missing"
rownames(pout)[is.na(rownames(pout))]="missing"
colnames(tout)[is.na(colnames(tout))]="missing"
colnames(pout)[is.na(colnames(pout))]="missing"

tout=as.data.frame(cbind(rownames(tout),as.data.frame.matrix(tout)))
names(tout)[1]=rowvar
if(!is.null(rowNames)) tout[,1]=rowNames

pout=as.data.frame(cbind(rownames(pout),as.data.frame.matrix(pout)))
names(pout)[1]=rowvar
if(!is.null(rowNames)) pout[,1]=rowNames

return(list(Counts=tout,Percent=pout))
}

##' Produces 2-way contingency tables, optionally with percentages, exports them to a spreadsheet, and saves the file.
##' 
##' This function produces two identically-sized tables side by side, one with the cross-tabulated counts of unique values of \code{rowvar, colvar}
##' and the other with percentages, calculated either by row (\code{sumby=1}, default) or column (\code{sumby=2}).
##' Row and column margins are also automatically produced. ##' Tables are automatically saved to the file associated with the \code{wb} spreadsheet object. 
##' 
##' There is an asymmetry between rows and columns, because the tables are converted to data frame in order for \code{\link{writeWorksheet}} to export them.

##' If you want to avoid the percentage tables on the side, choose \code{percents=FALSE}. If you also want no margins, just use the simpler function \code{\link{XLgeneric}}.
##' 

##' @note The worksheet \code{sheet} does not have to pre-exist; the function will create it if it doesn't already exist.
#' 
##' @note If \code{sheet} exists, it will be written into - rather than completely cleared and rewritten de novo. Only existing data in individual cells that are part of the exported tables' target range will be overwritten. If you do want to clear an existing sheet while exporting the new tables, set \code{purge=TRUE}. This behavior, and the usage of \code{purge}, are the same across all \code{table1xls} export functions.
##' 
##' 
##' @title Two-way Contingency Tables exported to a spreadsheet
##'
##' @param wb a \code{\link[XLConnect]{workbook-class}} object
##' @param sheet numeric or character: a worksheet name (character) or position (numeric) within \code{wb}.
##' @param rowvar vector: categorical variable (logical, numeric, character, factor, etc.) for the table's rows
##' @param colvar vector: categorical variable (logical, numeric, character factor, etc.) for the table's columns
##' @param sumby whether percentages should be calculated across rows (1, default) or columns (2).
##' @param rowTitle character: the title to be placed above the row name column (default empty string)
##' @param rowNames,colNames character vector of row and column names. Default behavior (\code{NULL}): automatically determined from data
##' @param ord numeric vector specifying row-index order in the produced table. Default (\code{NULL}) is no re-ordering.
##' @param row1,col1 numeric: the first row and column occupied by the table (title included if relevant).
##' @param title character: an optional overall title to the table. Default (\code{NULL}) is no title.
##' @param header logical: should a header row with the captions "Counts:" and "Percentages:" be added right above the tables? (default \code{FALSE}; ignored in any case if \code{percents=FALSE})
##' @param purge logical: should \code{sheet} be created anew, by first removing the previous copy if it exists? (default \code{FALSE})
##' @param digits numeric: how many digits (after the decimal point) to show in the percents table?
##' 
##' @param useNA How to handle missing values. Passed on to \code{\link{table}} (see help on that function for options).
##' @param percents logical: would you like only a count table (\code{FALSE}), or also a percents table side-by-side with the the count table (\code{TRUE}, default)?  
##' 
##' @return The function returns invisibly, after writing the data into \code{sheet}.
##' @example inst/examples/Ex2way.r 
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
##' @seealso Uses \code{\link{writeWorksheet}} to access the spreadsheet. See \code{\link{setStyleAction}} to control the output style. If interested in one-way tables, see \code{\link{XLoneWay}}.
##' @note This function uses the internal function \code{fancytab2} which produces 2-way tables with counts, percentages and margins. 

##' @export

XLtwoWay<-function(wb,sheet,rowvar,colvar,sumby=1,rowTitle="",rowNames=NULL,colNames=NULL,ord=NULL,row1=1,col1=1,title=NULL,header=FALSE,purge=FALSE,digits=1,useNA='ifany',percents=TRUE)
{

if(purge) removeSheet(wb,sheet)
if(!existsSheet(wb,sheet)) createSheet(wb,sheet)

### Producing counts and percents table via the internal function 'fancytab2'
tab=fancytab2(rowvar,colvar,sumby=sumby,rowvar=rowTitle,rowNames=rowNames,digits=digits,missings=useNA)


if(!is.null(title))  ### Adding a title
{
  XLaddText(wb,sheet,text=title,row1=row1,col1=col1)
  row1=row1+1
}


if (is.null(ord)) ord=1:dim(tab$Counts)[1]
if (!is.null(colNames)) 
{
    names(tab$Counts)[-1]=colNames
    names(tab$Percent)[-1]=colNames
}

widt=dim(tab$Counts)[2]+1
if(percents && header)  ### adding headers indicating 'counts' and 'percents'
{
  XLaddText(wb,sheet,"Counts:",row1=row1,col1=col1)
  XLaddText(wb,sheet,"Percent:",row1=row1,col1=col1+widt+1)
  row1=row1+1
}
writeWorksheet(wb,tab$Counts[ord,],sheet,startRow=row1,startCol=col1)
if(percents) writeWorksheet(wb,tab$Percent[ord,],sheet,startRow=row1,startCol=col1+widt+1)



setColumnWidth(wb, sheet = sheet, column = col1:(col1+2*widt+3), width=-1)
saveWorkbook(wb)

}  ### Function end


##' Univariate Statistics Exported to Excel
##' 
##' Calculates univariate summary statistics (optionally stratified), exports the formatted output to a spreadsheet, and saves the file.
##'
##' This function evaluates up to 2 univariate functions on the input vector \code{calcvar}, either as a single sample, or grouped by strata defined via \code{rowvar}. It produces a single-column table, each cell containing the formatted results from the two functions. The table is exported to a spreadsheet and the file is saved. 
##' 
##' The cell can be formatted to show a combined result, e.g. "Mean (SD)" which is the default. Tne function is quite mutable: both \code{fun1$fun, fun2$fun} and the strings separating their formatted output can be user-defined. The functions can return either a string (i.e., a formatted output) or a number that will be interpreted as a string in subsequent formatting.
##' The default calls \code{\link{roundmean},\link{roundSD}} and prints the summaries in \code{"mean(SD)"} format.
##' 
##' See the \code{\link{XLtwoWay}} help page, for behavior regarding new-sheet creation, overwriting, etc.

##' @return The function returns invisibly, after writing the data into \code{sheet} and saving the file.
##' 
##' @author Assaf P. Oron \code{<assaf.oron.at.seattlechildrens.org>}
##' @seealso Uses \code{\link{writeWorksheet}} to access the spreadsheet, \code{\link{rangeString}} for some utilities that can be used as \code{fun1$fun,fun2$fun}. For one-way (univariate) contingency tables, \code{\link{XLoneWay}}. 
##' 

##' 
##' @example inst/examples/ExUnivar.r 

##' @param wb a \code{\link[XLConnect]{workbook-class}} object
##' @param sheet numeric or character: a worksheet name (character) or position (numeric) within \code{wb}.
##' @param calcvar vector: variable to calculate the statistics for (usually numeric, can be logical).
##' @param rowvar vector: categorical variable to stratify \code{calcvar}'s summaries over. Default behavior if left unspecified, is to calculate overall summaries with the row title "All". 
##' @param fun1,fun2 two lists describing the utility functions that will calculate the statistics. Each list has a \code{fun} component for the function, and a \code{name} component for its name as it would appear in the column header.
##' @param seps character vector of length 3, specifying the formatted separators before the output of \code{fun1$fun}, between the two outputs, and after the output of \code{fun2$fun}. Default behavior encloses the second output in parentheses. See 'Examples'.
##' @param title character: an optional overall title to the table. Default (\code{NULL}) is no title.
##' @param rowTitle character: the title to be placed above the row name column (default empty string)
##' @param rowNames character vector of row names. Default behavior (\code{NULL}): automatically determined from data
##' @param ord numeric vector specifying row-index order (i.e., a re-ordering of \code{rowvar}'s levels) in the produced table. Default (\code{NULL}) is no re-ordering.
##' @param row1,col1 numeric: the first row and column occupied by the table (title included if relevant).
##' @param purge logical: should \code{sheet} be created anew, by first removing the previous copy if it exists? (default \code{FALSE})
##' @param ... parameters passed on to \code{fun1$fun,fun2$fun}
##'
##' @export

XLunivariate<-function(wb,sheet,calcvar,rowvar=rep("All",length(calcvar)),fun1=list(fun=roundmean,name="Mean"),fun2=list(fun=roundSD,name="SD"),seps=c('',' (',')'),title=NULL,rowTitle="",rowNames=NULL,ord=NULL,row1=1,col1=1,purge=FALSE,...)
{ 
if(purge) removeSheet(wb,sheet)  
if(!existsSheet(wb,sheet)) createSheet(wb,sheet)

num1=tapply(calcvar,rowvar,fun1$fun,...)
num2=tapply(calcvar,rowvar,fun2$fun,...)
if (is.null(ord)) ord=1:length(num1)
if (is.null(rowNames)) rowNames=names(num1)

outdat=data.frame(cbind(rowNames,paste(seps[1],num1,seps[2],num2,seps[3],sep='')))
names(outdat)=c(rowTitle,paste(seps[1],fun1$name,seps[2],fun2$name,seps[3],sep=''))

if(!is.null(title))  ### Adding a title
{
  XLaddText(wb,sheet,text=title,row1=row1,col1=col1)
  row1=row1+1
}
writeWorksheet(wb,outdat[ord,],sheet,startRow=row1,startCol=col1)

setColumnWidth(wb, sheet = sheet, column = col1:(col1+3), width=-1)
saveWorkbook(wb)
}
  

