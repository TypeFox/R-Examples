#' missing or zero pattern structure.
#' 
#' Analysis of the missing or the zero patterns structure of a data set.
#' 
#' Here, one pattern defines those observations that have the same structure
#' regarding their missingness or zeros. For all patterns a summary is
#' calculated.
#' 
#' @aliases missPatterns zeroPatterns
#' @param x a data frame or matrix.
#' @return \item{groups }{List of the different patterns and the observation
#' numbers for each pattern} \item{cn }{the names of the patterns coded as
#' vectors of 0-1's} \item{tabcomb}{the pattern structure - all combinations of
#' zeros or missings in the variables} \item{tabcombPlus}{the pattern structure
#' - all combinations of zeros or missings in the variables including the size
#' of those combinations/patterns, i.e. the number of observations that belongs
#' to each pattern.} \item{rsum}{the number of zeros or missing values in each
#' row of the data set}
#' @author Matthias Templ. The code is based on a previous version from Andreas
#' Alfons and Matthias Templ from package VIM
#' @seealso \code{\link[VIM]{aggr}}
#' @keywords multivariate
#' @export
#' @examples
#' 
#' data(expenditures)
#' ## set NA's artificial:
#' expenditures[expenditures < 300] <- NA
#' ## detect the NA structure:
#' missPatterns(expenditures)
#' 
missPatterns <- function(x){
	# identification of the missing pattern structure 
	# Matthias Templ, Oct 10, 2011
	if(is.null(dim(x))) stop("the data set has to be consist of at least two variables")
	
	w <- is.na(x)
	tmp <- ifelse(is.na(x), 1, 0)  # 'ifelse' does not omit 'dim' attribute
	tmpC <- apply(tmp, 1, paste, collapse=":")
	tab <- table(tmpC)
	tabcomb <- sapply(names(tab), 
			function(x) as.integer(unlist(strsplit(x, ":", fixed=TRUE))), 
			USE.NAMES=FALSE)
	tabcomb <- if(is.null(dim(tabcomb))) as.matrix(tabcomb) else t(tabcomb)
	tabcomb <- ifelse(tabcomb==0,TRUE,FALSE)
	cn <- names(tab)
	groups <- sapply(cn, function(y){
				(which(tmpC %in% y))
			})
	## Karels beiden MUSS-Variablen ;-):
	csum <- lapply(groups, length)
	amountComb <- cbind(data.frame(tabcomb), csum=as.numeric(csum))
	rsum <- apply(w, 1, sum)
	## TODO: N variable dazu, + 2. zeilenweise, spaltenweise
	list(groups=groups, cn=cn, tabcomb=tabcomb, tabcombPlus=amountComb, rsum=rsum)
}

#' @rdname missPatterns
#' @export
zeroPatterns <- function(x){
	# identification of the zero pattern structure 
	# Matthias Templ, Oct 10, 2011
	if(is.null(dim(x))) stop("the data set has to be consist of at least two variables")
	
	w <- x == 0
	tmp <- ifelse(x==0, 1, 0)  # 'ifelse' does not omit 'dim' attribute
	tmpC <- apply(tmp, 1, paste, collapse=":")
	tab <- table(tmpC)
	tabcomb <- sapply(names(tab), 
			function(x) as.integer(unlist(strsplit(x, ":", fixed=TRUE))), 
			USE.NAMES=FALSE)
	tabcomb <- if(is.null(dim(tabcomb))) as.matrix(tabcomb) else t(tabcomb)
	tabcomb <- ifelse(tabcomb==0,TRUE,FALSE)
	cn <- names(tab)
	groups <- sapply(cn, function(y){
				(which(tmpC %in% y))
			})
	## Karels beiden MUSS-Variablen ;-):
	csum <- lapply(groups, length)
	amountComb <- cbind(data.frame(tabcomb), csum=as.numeric(csum))
	rsum <- apply(w, 1, sum)
	## TODO: N variable dazu, + 2. zeilenweise, spaltenweise
	list(groups=groups, cn=cn, tabcomb=tabcomb, tabcombPlus=amountComb, rsum=rsum)
}
