io2 <- function(...,data=NULL) {
# ... can be a sequence of lists of numeric sample vectors 
# with at least 2 sample vectors per list.
# Each list corresponds to a block level, and the samples
# within a block correspond to different (treatment) 
# group levels, although the same group levels (treatments) 
# may be used within different blocks.
# Or ... can be a list of such lists.
# Or ... can be a formula that specifies a
# a response (e.g., y), grouped by a treatment factor, e.g., g, 
# with a blocking factor, e.g., b, via formula y ~ g | b. 
# The variables y, g, b may be names of columns in a supplied
# data frame dat via data = dat (default NULL).
# When data = NULL is used the variables y, g, b should exist 
# in the calling environment.
# This breaks down y into blocks of responses via the block
# levels and within each block the responses are broken down into 
# different samples corresponding to the different levels of g.
# Or ... can be a data.frame with three column y, g, and b, 
# as described above.
xlist <- list(...)
cl <- match.call()  # gets a copy of the current call
if(is(xlist[[1L]], "formula")) {
	clstr <- unlist (strsplit(as.character(eval(cl[[2]])),NULL))
	if ("~" %in% clstr & "|" %in% clstr) {
			
			Y <- 	eval(as.name(strsplit(
				as.character(eval(cl[[2]]))," ")[[2]]),
				envir = data)
			G <- eval(as.name(strsplit(
				as.character(eval(cl[[2]]))," ")[[3]][1]),
				envir = data)
			B <- eval(as.name(strsplit(
				as.character(eval(cl[[2]]))," ")[[3]][3]),
				envir = data)
			data.sets <- listmake(
						Y,as.factor(G),as.factor(B))
	}
} else {	if(length(xlist) == 1 && is.list(...[[1]])){	
			data.sets <- list(...)[[1]]
		} else {
			if( length(list(...)) > 1 && 
				is.list(list(...)[[1]])){
				data.sets <- list(...)
			} 
		}
 	}
b <- length(data.sets)
if(b < 2) stop("less than 2 blocks\n")
num <- NULL
for( i in 1:b ){
	num <- c(num,unlist(lapply(data.sets[[i]],FUN=is.numeric)))
        if(length(data.sets[[i]]) < 2) stop("at least one block with less than 2 samples\n")
}
if( all(num) == FALSE ) stop("improper data in ... \n")
data.sets
}

