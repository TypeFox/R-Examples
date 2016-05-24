### return a data.frame giving brief summaries of data sets in packages
#  package:    a character vector giving the package(s) to look in
#  allClass:   include all classes of the item (TRUE) or just the last class (FALSE)
#  incPackage: include package name in result?
#  maxTitle: maximum length of data set Title

datasets <- function(package, allClass=FALSE, 
		incPackage=length(package) > 1,
		maxTitle=NULL) 
{
	# make sure requested packages are available and loaded
	pkgs <- .packages()
	for (i in seq_along(package)) {
		if (! package[i] %in% pkgs) 
			if (require(package[i], character.only=TRUE, quietly=TRUE))
				cat(paste("Loading package:", package[i], "\n"))
			else stop(paste("Package", package[i], "is not available"))
	}
	dsitems <- data(package=package)$results
	wanted <- if (incPackage) c('Package', 'Item','Title') else c('Item','Title')
	ds <- as.data.frame(dsitems[,wanted], stringsAsFactors=FALSE)
	# fix items with " (...)" in names, e.g., "BJsales.lead (BJsales)" in datasets
	ds$Item <- gsub(" .*", "", ds$Item)
	
	getDim <- function(x) {
		if (is.null(dim(get(x)))) length(get(x)) else paste(dim(get(x)), collapse='x')
	}
	getClass <- function(x) {
		cl <- class(get(x))
		if (length(cl)>1 && !allClass) cl[length(cl)] else cl
	}
	ds$dim <- unlist(lapply(ds$Item, getDim ))
	ds$class <- unlist(lapply(ds$Item, getClass ))
	if (!is.null(maxTitle)) ds$Title <- substr(ds$Title, 1, maxTitle)
	if (incPackage)
		ds[c('Package', 'Item','class','dim','Title')]
	else
		ds[c('Item','class','dim','Title')]
}
