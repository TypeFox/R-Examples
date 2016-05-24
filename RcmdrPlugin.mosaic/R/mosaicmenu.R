mosaicMenuItem <- function() {
	if (weHaveMosaicableTables()) {
		mosaicSelectDialog() 
	} else {
		Message("No tables available for this type of plot","warning")
	}
}		


weHaveMosaicableTables <- function() {
	length(Filter(function(x)(inherits(get(x),"table")|inherits(get(x),"ftable")|inherits(get(x),"structable")),ls(envir=globalenv())))>0
}
