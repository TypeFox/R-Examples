# TODO: Add comment
# 
# Author: ecor
###############################################################################

#
#> fruits <- c("one apple", "two pears", "three bananas")
#> str_replace(fruits, "[aeiou]", "-")
#[1] "-ne apple"     "tw- pears"     "thr-e bananas"
#> fruits
#[1] "one apple"     "two pears"     "three bananas"
#> fruits <- c( "apples and oranges and pears and bananas", "pineapples and mangos and guavas"
#				+ )
#> fruits
#[1] "apples and oranges and pears and bananas" "pineapples and mangos and guavas"        
#> str_split(fruits, " and ")
#[[1]]
#[1] "apples"  "oranges" "pears"   "bananas"
#
#[[2]]
#[1] "pineapples" "mangos"     "guavas"    
#
#> str_split(fruits, "!")
#[[1]]
#[1] "apples and oranges and pears and bananas"
#
#[[2]]
#[1] "pineapples and mangos and guavas"
#
#> fruits <- c( "!apples and oranges and ! pears and bananas", "pineapples and mangos and guavas"
#				+ )
#> str_split(fruits, "!")
#[[1]]
#[1] ""                        "apples and oranges and " " pears and bananas"     
#
NULL
#'
#' Collects all keywords contained in the 'getop.inpts' configuration files and their values in a data frame object. 
#' 
#' @param wpath working directory containing GEOtop files
#' @param inpts.file name of the GEOtop configuration file. Default is \code{"geotop.inpts"}
#' @param comment comment indicator charcater. Default is \code{"!"}
#' @param no.comment string indicatos read as comment ones by GEOtop but they do not indicate comments by "geotopbricks" package. 
#' @param exceptions string vector. If keywords contain an element of this vector, the blank spaces in Value \code{" "} will not be removed. 
#' @param warn logical argument of \code{\link{readLines}}. Default is \code{FALSE}.
#' @param ... further arguments of \code{\link{readLines}} 
#' 
#' @export 
#' 
#' @return a data frame with two columns: \code{Keyword} and \code{Value}
#' 
#' @import stringr 

#' @seealso \code{\link{get.geotop.inpts.keyword.value}}
#' 

declared.geotop.inpts.keywords <- function(wpath,inpts.file="geotop.inpts",comment="!",exceptions="Date",warn=FALSE,no.comment=c("!>!","!>>!"),...) {

	
	if (!is.null(wpath)) {
					
		file <- paste(wpath,inpts.file,sep="/") 
		
	} else {
		
		file <- inpts.file
		
	}
	
	
#	if (str_sub(file,1,3)=='ssh' | str_sub(file,1,5)=='plink') {
#		
#		file <- pipe(file) # added line according to http://stackoverflow.com/posts/2226880/edit
#		open <- TRUE
#	}	else {
#		print('qui')
#		file <- file(file)
#		open <- FALSE
#	} ## commented line (to be removed) by ec 2014-05-20
#	file <- file(file)
	x <- readLines(file,warn=warn,...)
	
#	if (open) close(file)
	
	l <- nchar(x)
	
	x <- x[l>2] 
	
	if (length(no.comment)>0) {
		
		xout <- x 
		for (itnc in no.comment) {
			
		##print(xout[35:38])
		xout <- unlist(lapply(X=xout,FUN=str_replace, pattern=itnc, replacement=""))
		##print(xout[35:58])
		}
	} else {
		
		xout <- x
	}
	xout <- str_split(xout,comment)

	
	for (i in 1:length(xout)) {
		
		x[i] <- xout[[i]][1]
		
	}
	l <- nchar(x)
	x <- x[l>2] 
	x <- str_replace_all(x,c("\t"),"")
	x <- x[!is.na(x)]
	
	
#	x <- str_replace_all(x,c(" "),"")
	split="="
	x <- x[str_detect(x,split)]
	
	xout <- str_split(x,split)
	
	out <- as.data.frame(array(NA,c((length(xout)),2)))
	
	names(out) <- c("Keyword","Value")
	for (i in 1:length(xout)) {
		
		
		out$Keyword[i] <- xout[[i]][1] 
		out$Value[i] <- xout[[i]][2]
	
		out$Keyword <- str_replace_all(out$Keyword,c(" "),"")		
	}	
	
	out$Keyword <- str_replace_all(out$Keyword,c(" "),"")	
	
	exc <- array(FALSE,length(out$Keyword))
	
	
	for (i in 1:length(exceptions)) {
		
		exc <- str_detect(out$Keyword,exceptions[i]) | exc
		
		
	}
	
	out$Value[!exc] <- str_replace_all(out$Value[!exc],c(" "),"")
	
	
	out <- out[!is.na(out$Keyword),] ### correction by ec on 20131104 
	out <- out[out$Keyword!="",] 
	out <- out[out$Value!="",] 
	
	
	return(out)
	
}


