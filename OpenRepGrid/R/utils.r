# Variance, standard deviation and covariance without Bessel's correction

covpop <- function(x, y, na.rm=TRUE){
	x <- unlist(x)
	y <- unlist(y)
  if (na.rm) {                    # delete missings
    index <- is.na(x | is.na(y))
    x <- x[!index]
    y <- y[!index]
  }
	n <- length(x)
	((n-1)/n) * cov(x=x, y=y)			# undo Bessel's correction
}

varpop <- function(x, na.rm=FALSE){
  covpop(x=x, y=x, na.rm=na.rm)			# undo Bessel's correction
}

sdpop <- function(...){
	sqrt(varpop(...))
}


# factorial function 
# wrapper for convenience
fac <- function (x) gamma(1 + x)

joinString <- function(x) 
  paste(unlist(x), sep="", collapse=" ")
  
trimBlanksInString <- function(x) 
  sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", x, perl=TRUE)


baseSplitStringInt <- function(text, availwidth=1, cex=1)			# function to split text in base graphics
{				
	if (is.expression(text)){			# expressions cannot be split
		return(text)
		break
	}
	if (identical(text, NULL)) text <- ""
	if (identical(text, NA)) text <- ""
	if (identical(text, character(0))) text <- ""
	if (text == ""){
		return(paste(text))
		break
	}

	strings <- strsplit(as.character(text), " ")[[1]]
	if (length(strings) == 1){
		return(paste(strings))
		break
	}
	newstring <- strings[1]
	linewidth <- strwidth(newstring, cex = cex)
	gapwidth <-  strwidth(" ", cex = cex) 

	for (i in 2:length(strings)) {
		width <- strwidth(strings[i], cex = cex)
		if (linewidth + gapwidth + width < availwidth){
			sep <- " "
			linewidth <- linewidth + gapwidth + width
		} else {
			sep <- "\n"
			linewidth <- width
		}
		newstring <- paste(newstring, strings[i], sep=sep)
	}
	newstring
}


baseSplitString <- function(text, availwidth=1, cex=1){
  as.vector(sapply(text, baseSplitStringInt, 
                   availwidth=availwidth, cex=cex))
}


# makeStandardRangeColorRamp() creates color ramp for supplied colors that takes 
# values between [0,1] and returns a hex color value
#
makeStandardRangeColorRamp <- function(colors, na.col="#FFFFFF", ...){
	ramp <- colorRamp(colors, ...)
	function(x){
    is.na(x) <- is.na(x)      # convert NaN values to NA
    na.index <- is.na(x)
    x[na.index] <- 0          # overwrite so color can be determined   
    x <- ramp(x)              # actual color calculation
    col <- rgb(x[, 1], x[, 2], x[, 3], maxColorValue = 255)	  
    col[na.index] <- na.col   # replace na indices with default NA color
    col                   
	}	
}


#' modifyListNull
#' 
#' TODO: a modified version of modifyList that does not overwrite elements
#' if they are NULL in the supplied list
#'
#' @param   x
#' @param   val
#' @return  list
#' @noRd
#'
modifyListNull <- function (x, val) 
{
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) 
            Recall(x[[v]], val[[v]])
			else if(!is.null(val[[v]])){					# this part was extended to check if element is NULL
				val[[v]]
			} else x[[v]]
    }
    x
}
#l1 <- list(a=1, b=2)
#l2 <- list(a=NULL, b=3)
#modifyListNull(l1, l2)
#modifyList(l1,l2)
#modifyListNull(l2, l1)

#' modifyListNA
#' 
#' TODO: a modified version of modifyList that does not overwrite elements 
#' if they are NA in the supplied list
#' @param   x
#' @param   val
#' @return  list
#' @noRd
#'
modifyListNA <- function (x, val) {
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    for (v in names(val)) {
        x[[v]] <- if (v %in% xnames && is.list(x[[v]]) && is.list(val[[v]])) {
				Recall(x[[v]], val[[v]])
			} else if (!is.na(val[[v]])) {					# this part was extended to check if element is NULL
				val[[v]]
			} else x[[v]]
    }
    x
}
#l1 <- list(a=1, b=2)
#l2 <- list(a=NA, b=3)
#modifyListNA(l1, l2)
#modifyList(l1,l2)
#modifyListNA(l2, l1)

#l1 <- list(t=list(a=1, b=2))
#l2 <- list(t=list(a=NA, b=3))
#modifyListNA(l1, l2)
#modifyList(l1,l2)
#modifyListNA(l2, l1)


###############################################################################
#' bring vector values into ring form
#'
#' the values of a vector that are outside of a given range are remapped
#' to the values of the range. This function is useful for loops over rows and
#' columns of a matrix if the 
#' 
#' @param x       vector
#' @param upper   upper limit of range (lower is one. TODO: maybe extend???)
#' @return vector
#' @export
#' @keywords internal
#' @examples \dontrun{
#'    ring(1:10, 3)
#'
#'    m <- matrix(1:12, ncol=4)
#'    for(i in 1:12)
#'      print(m[ring(i, 3), map(i, 4)])
#' }
#'
ring <-  function(x, upper){
  res <- x %% upper
  res[res == 0] <- upper
  res
}


###############################################################################
#' map a value onto others
#' 
#' @param x      vector
#' @param each   number of cuts
#' @return vector
#' @export
#' @keywords internal
#' @examples \dontrun{
#'    map(1:10, 3)
#'
#'    m <- matrix(1:12, ncol=4)
#'    for(i in 1:12)
#'      print(m[ring(i, 3), map(i, 4)])
#' }
#'
map <- function(x, each){
  ceiling(x/each)
}


###############################################################################
#' order one vector by another
#'
#' small wrapper to order one vector by another, hardly worth a function
#' 
#' @param x   vector
#' @param y   vector
#' @return    vector
#' @export
#' @keywords internal
orderBy <- function(x,y) y[order(x)]


###############################################################################
#' make ascending and descending vector
#'
#' along a given length n make ascending indices until reaching
#' the midpoint and descend afterwards again.
#' 
#' @param n       \code{integer} The length of the indexes
#' @param type    (integer, default=1). If 1 the cascade index is returned. 
#'                2 returns the index of left and right side, 3 returns the length
#'                of the left and right index vector
#' @return        vector (type 1 or 3) or list (type 2)
#' @export
#' @keywords internal
#' @examples \dontrun{
#'    for(n in 1:10) 
#'        print(cascade(n)) 
#' }
cascade <- function(n, type=1){
  if (type == 2){
    list( left=(1:n)[0:floor(n/2)], 
          right=rev((n:1)[0:ceiling(n/2)]) )
  } else if (type == 3){
    c( left=length((1:n)[0:floor(n/2)]), 
       right=length((n:1)[0:ceiling(n/2)]) )
  } else {
    c((1:n)[0:floor(n/2)], rev((1:n)[0:ceiling(n/2)]))
  }
}


# insertAt kreiert die Indizes für das ineinanderfügen von zwei Vektoren, Listen etc.
# index.base		Index des Objekts in das eingefügt werden soll (meist 1,2,3 etc.)
# index.insert		Index der Stellen an denen ein Objekt eingefügt werden soll
#
# 1 2 3 4		1  2  3  4
#1   3		   1           5
#   2 3 4 5
# 1    4
#   2 3 5 6
# 1    4

#' insertAt
#' 
#' TODO: a modified version of modifyList that does not overwrite elements 
#' if they are NA in the supplied list
#' @param   x
#' @param   val
#' @return  list
#' @noRd
insertAt <- function(index.base, index.insert, side="pre"){
	if(!side %in% c("pre", "post")) 									# Integrity Checks
		stop("side must be a a string. It can take the values 'pre' or 'post'")					
	res <- list(index.base=index.base, index.insert=index.insert)
	for(i in seq_along(index.insert)){
		at <- index.insert[i]
		if(side=="pre"){												# VOR der benannten Position at einfügen
			index.base <- index.base + (index.base >= at)				# Alle Indizes größer-gleich at werden um eines erhöht
			options(warn=-1)											# in case index.base=numeric(0) warnings gets generated at max()
			index.insert <- index.insert + ((index.insert > at) & 
								any(index.insert[seq_along(index.insert) > i] <= max(index.base)))					
			options(warn=0)
		}	
		if(side=="post"){												# NACH der benannten Position at einfügen
			index.base <- index.base + (index.base > at)				# Alle Indizes größer als at werden um eines erhöht
			options(warn=-1)											# in case index.base=numeric(0) warnings gets generated at max()
			index.insert <- index.insert + ((index.insert >= at) & 
								any(index.insert[seq_along(index.insert) >= i] <= max(index.base)))					
			options(warn=0)
		}
	}	
	c(res, list(index.base.new=index.base, index.insert.new=index.insert))
}



# insertAt(1:4, c(1,3))
# insertAt(c(1,2,3,4), c(1,3), side="pre")
# insertAt(c(1,2,3,4), c(1,2), side="pre")
# insertAt(c(1,2,3,4), c(4,5), side="pre")
# insertAt(c(1,2,3,4), c(5,6), side="pre")
# insertAt(c(1,2,3,4), c(1,2,3,5,6), side="pre")
# insertAt(1:4, 5:8)
# insertAt(numeric(0), 1:2)
# insertAt(numeric(0), c(1,3))

# insertAt(c(1,2,3,4), c(1,3), side="post")
# insertAt(c(1,2,3,4), c(1,2), side="post")
# insertAt(c(1,2,3,4), c(4,5), side="post")
# insertAt(c(1,2,3,4), c(5,6), side="post")
# insertAt(c(1,2,3,4), c(1,2,3,5,6), side="post")
# insertAt(1:4, 5:8, side="post")
# insertAt(numeric(0), 1:2, side="post")
# insertAt(numeric(0), c(1,3), side="post")



# There was once question on r-help asking if apply could be used with a progress bar.
# The plyr package provides several apply like functions also including progress bars,
# so one could have a look here and use a plyr function instead of apply if possible.
# Anyway, here comes a wrapper for apply/lapply that has a progressbar.

# Here is a wrapper for a function passed to apply that will create a text progress bar

# STATUS: WOKRING, but only tested once or twice, tested with ?apply examples
# ISSUES/TODO: MARGIN argument not always correct when vector like c(1,2) is used





#' apply with a progress bar
#'
#' Can be used like standard base:::apply. The only thing 
#' it does is create an additional progress bar.
#'
#' @param X         see ?apply for parameter explanation
#' @param MARGIN    see ?apply
#' @param FUN       see ?apply
#' @param ...       see ?apply
#' @return          see ?apply
#' @seealso         \code{\link{apply}}
#' @author          Mark Heckmann
#' @export
#' @keywords        internal
#' @examples \dontrun{
#'
#'    apply_pb(anscombe, 2, sd, na.rm=TRUE)
#'
#'    # larger dataset
#'    df <- data.frame(rnorm(30000), rnorm(30000))
#'    head(apply_pb(df, 1, sd))
#'
#'    # performance comparison
#'    df <- data.frame(rnorm(90000), rnorm(90000))
#'    system.time(apply(df, 1, sd))
#'    system.time(apply_pb(df, 1, sd))
#'
#' }
#'
apply_pb <- function(X, MARGIN, FUN, ...)
{
	env <- environment()                                      # this environment
	pb_Total <- sum(dim(X)[MARGIN])                           # get mex value for progress bar
	counter <- 0                                              # make counter variable
	pb <-  txtProgressBar(min = 0, max = pb_Total, style = 3) # make progress bar
		
	# wrapper around FUN
	wrapper <- function(...){
		curVal <- get("counter", envir =  env)              # get counter value 
		assign("counter", curVal +1 ,envir= env)            # and increment it by one
		setTxtProgressBar(get("pb", envir= env), curVal +1) # update progress bar
		FUN(...)
	}
	res <- apply(X, MARGIN, wrapper, ...)  # apply wrapper with apply  
	close(pb)                              # close progress bar
	res
}

# apply_pb(anscombe, 2, sd, na.rm=TRUE)

# large dataset
# df <- data.frame(rnorm(30000), rnorm(30000))
# head(apply_pb(df, 1, sd))



#' lapply with a progress bar
#'
#' Can be used like standard base:::lapply. The only thing 
#' it does is create an additional progress bar.
#'
#' @param X           see ?lapply for parameter explanation
#' @param FUN         see ?lapply 
#' @param ...         see ?lapply 
#' @return list       see ?lapply
#' @seealso           \code{\link{lapply}}
#' @author            Mark Heckmann
#' @export
#' @keywords          internal
#' @examples \dontrun{
#'
#'    l <- sapply(1:20000, function(x) list(rnorm(1000)))
#'    lapply_pb(l, mean)
#'
#' }
#'
lapply_pb <- function(X, FUN, ...)
{
	env <- environment()                                      # this environment
	pb_Total <- length(X)                                     # get max value for progress bar
	counter <- 0                                              # make counter variable
	pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)  # make progress bar
	
	# wrapper around FUN
	wrapper <- function(...){
		curVal <- get("counter", envir = env)              	# get counter value 
		assign("counter", curVal +1 ,envir=env)            	# and increment it by one
		setTxtProgressBar(get("pb", envir=env), curVal +1)  # update progress bar
		FUN(...)
	}
	res <- lapply(X, wrapper, ...)  # use wrapper with lapply
	close(pb)                       # close progress bar
	res                                	  
}

# l <- lapply(1:20000, function(x) list(rnorm(1000)))
# head(lapply_pb(l, mean))



#' sapply with a progress bar
#'
#' Can be used like standard base:::sapply. The ionly thing 
#' it does is create an additional progress bar.
#'
#' @param X           see ?sapply for parameter explanation
#' @param FUN         see ?sapply 
#' @param ...         see ?sapply 
#' @return list       see ?sapply
#' @seealso           \code{\link{sapply}}
#' @author            Mark Heckmann
#' @export
#' @keywords          internal
#' @examples \dontrun{
#'
#'    l <- sapply(1:20000, function(x) list(rnorm(1000)))
#'    head(sapply_pb(l, mean))
#'
#'    # performance comparison
#'    l <- sapply(1:20000, function(x) list(rnorm(1000)))
#'    system.time(sapply(l, mean))
#'    system.time(sapply_pb(l, mean))
#'
#' }
sapply_pb <- function(X, FUN, ...)
{
	env <- environment()                                      # this environment
	pb_Total <- length(X)                                     # get max value for progress bar
	counter <- 0                                              # make counter variable
	pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)  # make progress bar
		
	# wrapper around FUN
	wrapper <- function(...){
		curVal <- get("counter", envir = env)              	# get counter value 
		assign("counter", curVal +1 ,envir=env)            	# and increment it by one
		setTxtProgressBar(get("pb", envir=env), curVal +1)  # update progress bar
		FUN(...)
	}
	res <- sapply(X, wrapper, ...)   	# use wrapper with sapply
	close(pb)                       	# close progress bar
	res
}


#' reverse a string
#' 
#' reverses the strings of a vector, i.e. c("ABC", "abc")
#' becomes c("CBA", "cba")
#'
#' @param x   a string or a vector of strings
#' @return vector  a string or vector of strings with reversed chars
#' @references
#'    \url{http://www.mail-archive.com/r-help@@r-project.org/msg102759.html}
#' @export
#' @keywords internal
#' @examples
#' strReverse(c("ABC", "abc"))
strReverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), 
         paste, collapse="")
}


#' trim vector to lower or upper value 
#'
#' the range a value may take is resticted to by an upper and 
#' lower boundary. If it excedd the boundary the value is replaced
#' by the boundary value or alternatively by NA
#'
#' @param x         numeric vector
#' @param minmax    minimal and maximal possible value (default c(-Inf, Inf) 
#'                  i.e. no trimming occurs)
#' @param na        Use NAs for replacing values that are out of range
#' @return vector   vector whose elemenets that are out of range are replaced
#' @export
#' @keywords internal
#' @examples
#' trim_val(30)
#' trim_val(30, c(10,20))
#'
trim_val <- function(x, minmax=c(-Inf, Inf), na=FALSE){
  if(na){
    x[x < minmax[1]] <- NA
    x[x > minmax[2]] <- NA
  }
  else {
    x[x < minmax[1]] <- minmax[1]
    x[x > minmax[2]] <- minmax[2] 
  }
  x
}




#' recycle vector to given length
#'
#' @param vec       vector to be recycled
#' @param length    integer or vector. integer determines length of 
#'                  recycling. If a vector is provided the length of the 
#'                  vector is used.
#' @param na.fill   Use NAs for filling up to given length (default=FALSE)
#' @return vector
#' @note If 2nd argument is a vector, the first argument is recycled
#' to the length of the second vector. Instead of recycling the vector can 
#' also be added extra NAs if the length argument is smaller than the 
#' number of elements from vec, vec is cut off to make it usable for 
#' many purposes.
#' 
#' @export
#' @keywords internal
#' @examples
#' recycle(c(1,2,3), 7)
#' recycle(c(1,2,3), letters[1:7])
#' recycle(c(1,2,3), 7, na.fill=TRUE)
#' recycle(1, letters[1:3], na.fill=TRUE)
#' recycle(letters[1:3], 7)
#' recycle(letters[1:3],  letters[1:7])
#' recycle(letters[1:40],  letters[1:7])		# vec is cut off
recycle <- function(vec, length, na.fill=FALSE){
	if (!is.vector(vec) & !is.vector(length))
		stop("vec and length must be vectors. length may also be an integer")
	if (!is.numeric(length) & is.vector(length))	# both vectors
		length <- length(length)
	if (is.vector(length) & length(length) > 1L)	# is length a vector longer than 1
		length <- length(length)							      # then get length of vector
	if (!na.fill) {
		newvec <- rep(vec, ceiling(length / length(vec)))	  # enlarge vector by recycling
	} else {
		newvec <- c(vec, rep(NA, length * 
		            (ceiling(length / length(vec)) - 1L)))  # fill up with NAs
	}
	newvec[1L:length]
}



#' variation of recycle that recycles one vector x or y to the length of the 
#' longer one
#'
#'
#' @param x         vector to be recycled if shorter than y
#' @param y         vector to be recycled if shorter than x
#' @param na.fill   Use NAs for filling up to given length (default=FALSE)
#' @return list     a list containing the recycled x vector as first and 
#'                  the recycled y vector as second element             
#' @export
#' @keywords internal
#' @examples
#' recycle2(1:10, 1:3)
#' recycle2(1, 1:5)
#' recycle2(1, 1:5, na.fill=TRUE)
#' recycle2(1:5, 5:1)    # vectors unchanged
recycle2 <- function(x, y, na.fill=FALSE){
  len.x <- length(x)
  len.y <- length(y)
  if(len.x < len.y)
    x <- recycle(x, len.y, na.fill) 
  else if (len.x > len.y)
    y <- recycle(y, len.x, na.fill)
  list(x=x,y=y)
}




#' generate a random words
#' 
#' randomWords generates a vector of random words taken from a small 
#' set of words
#' @param n number of words to be generated (integer)
#' @return a string with n words (if length is not constrained)
#' @export
#' @keywords internal
#' @examples
#' randomWords(10)  # 10 random words
randomWords <- function(n){
	if (! is.numeric(n))
	  stop("n must be an integer")
	words <- c( "the", "novel", "depicts", "Harry", "as", "an", "essentially",
	            "good", "man", "who", "is", "forced", "into", "blackmarket", 
	            "activity", "by",	"economic", "forces", "beyond", "his", 
	            "control", "initially", "his", "fishing", "charter", 
	            "customer", "Mr.", "Johnson", "tricks", "Mark", "by", 
	            "slipping", "away", "without", "paying", "any", "of", "the",
	            "money", "he", "owes", "him", "Brownstone", "then", "flees", 
	            "back", "to", "the", "mainland", "by", "airplane", "before", 
	            "he", "realizes", "what", "has", "happened", "I", "she")
	sample(words, n, replace=TRUE)
}


#' generate a random sentence with n words
#'
#' @param n   number of word in sentence
#' @param maxchar   maximal number of characters per sentence. Note that whole 
#'                  words (not part of words) are excluded if the maximal number 
#'                   is exceeded.
#' @return a string with n words (if length is not constrained)
#' @export
#' @keywords internal
#' @examples  
#' randomSentence(10)   # one random sentence with 10 words
randomSentence <- function(n, maxchar=Inf){
	x <- paste(randomWords(n), collapse=" ")
	x.split <- strsplit(x, " ")[[1]]
	chars <- as.vector(sapply(x.split, nchar))
	paste(unlist(x.split[cumsum(chars) < maxchar]), collapse = " ")
}


#' generate n random sentences with a given or random number of words
#'
#' @param n         number of sentences to be generate (integer)
#' @param nwords    number of words per sentence. If vector each sentence
#'           lengths is randomly drawn from the vector
#' @param maxchar   maximal number of characters per sentence. Note that whole 
#'           words (not part of words) are excluded if the maximal number 
#'          is exceeded.
#' @return a vector with n random sentences
#' @export
#' @keywords internal
#' @examples
#' randomSentences(5, 10)     # five random sentences with ten words each
#' randomSentences(5, 2:10)   # five random sentences between two and ten words
randomSentences <- function(n, nwords, maxchar=Inf){
  sapply(sample(nwords, n, replace = TRUE), 
          randomSentence, maxchar = maxchar)
}


#' find the order of a string vector so it will match the order of another
#'
#' @param x   a vector of strings
#' @param y   a vector of strings
#' @return  a vector of strings
#' @export
#' @keywords internal
#' @examples \dontrun{
#'    a <- c("c", "a", "b")
#'    b <- c("b", "c", "a")
#'    index <- orderByString(a, b)    # to order b like a needs what indexes?
#'    index
#'    b[index]
#' }
#'
orderByString <- function(x, y){
  if (!all(x %in% y))
    stop("vector x and y do not contain the same (differently ordered) elements")
  index <- order(order(x))      # reconversion index from sorted to old order 
  order(y)[index]
}



### Thanks to Marc Schwartz for supplying the code for the Somer's d measure

# Calculate Concordant Pairs in a table
# cycle through x[r, c] and multiply by
# sum(x elements below and to the right of x[r, c])
# x = table
concordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND > c)
  # for each matrix[r, c]
  mat.lr <- function(r, c)
  { 
    lr <- x[(r.x > r) & (c.x > c)]
    sum(lr)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.lr, r = r.x, c = c.x))
}

# Calculate DIScordant Pairs in a table
# cycle through x[r, c] and multiply by
# sum(x elements below and to the left of x[r, c])
# x = table
discordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND < c)
  # for each matrix[r, c]
  mat.ll <- function(r, c)
  { 
    ll <- x[(r.x > r) & (c.x < c)]
    sum(ll)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.ll, r = r.x, c = c.x))
}


# Calculate Somers' d
# Return 3 values:
# 1. Sd C~R
# 2. Sd R~C
# 3. Sd Symmetric (Mean of above)
# x = table
calc.Sd <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)
  n <- sum(x)
  SumR <- rowSums(x)
  SumC <- colSums(x)

  Sd.CR <- (2 * (c - d)) / ((n ^ 2) - (sum(SumR ^ 2)))
  Sd.RC <- (2 * (c - d)) / ((n ^ 2) - (sum(SumC ^ 2)))
  Sd.S <- (2 * (c - d)) / ((n ^ 2) - (((sum(SumR ^ 2)) + (sum(SumC ^ 2))) / 2))

  Sdlist <- list(Sd.CR, Sd.RC, Sd.S)
  names(Sdlist) <- c("Sd.CR", "Sd.RC", "Sd.S")

  Sdlist
}

## example from Kaehler book, p.123 table, p.129 results
# m <- matrix(c(4,6,0,11,146,22,2,20,39), 3)
# calc.Sd(m)    # correct



# ellipse and angle code from: Dr P.D.M. Macdonald
# http://www.math.mcmaster.ca/peter/s4c03/s4c03_0506/classnotes/DrawingEllipsesinR.pdf

# draw an ellipse
#
# 
ellipse <- function (hlaxa = 1, hlaxb = 1, theta = 0, xc = 0, yc = 0, 
                    newplot = F, npoints = 100, ...)
{
  a <- seq(0, 2 * pi, length = npoints + 1)
  x <- hlaxa * cos(a)
  y <- hlaxb * sin(a)
  alpha <- angle(x, y)
  rad <- sqrt(x^2 + y^2)
  xp <- rad * cos(alpha + theta) + xc
  yp <- rad * sin(alpha + theta) + yc
  if (newplot)
    plot(xp, yp, type = "l", ...)
  else lines(xp, yp, ...)
  invisible()
}

angle <- function (x, y)
{
  angle2 <- function(xy) {
    x <- xy[1]
    y <- xy[2]
    if (x > 0) {
      atan(y/x)
    } else {
      if (x < 0 & y != 0) {
        atan(y/x) + sign(y) * pi
      }
      else {
        if (x < 0 & y == 0) {
          pi
        }
        else {
          if (y != 0) {
            (sign(y) * pi)/2
          }
          else {
            NA
          }
        }
      }
    }
  }
  apply(cbind(x, y), 1, angle2)
}


###############################################################################
###                           FORMATTING                                   ####
###############################################################################

#' Format a matrix and add index column.
#'
#' @param x         A matrix onbject.
#' @param rnames    Row names.
#' @param cnames    Column names.
#' @param pre.index Whether to make index prefix for rows and column names.
#' @param indexcol  Whether to make an index column.
#' @param diag      Whether to show diagonal.
#' @param mode      Whether to show upper (mode=1), lower (mode=2) 
#'                  or both triangles (mode=0) of the matrix.
#'
#' @author          Mark Heckmann
#' @keywords        internal
#' @export
#'
formatMatrix <- function(x, rnames=rownames(x), pre.index=c(T,F),
                         cnames=seq_len(ncol(x)), indexcol=F, digits=2,  
                         diag=F, mode=1)
{
  blanks <- paste(rep(" ", digits + 2), collapse="")
  if (mode == 1)
    x[lower.tri(x, diag=!diag)] <- blanks
  if (mode == 2)
    x[upper.tri(x, diag=!diag)] <- blanks
  
  if (pre.index[1])
    rnames <- paste(seq_len(nrow(x)), rnames) 
  if (pre.index[2])
    cnames <- paste(seq_len(ncol(x)), cnames)
  if (indexcol) {
    rownames(x) <- rnames
    x <- addIndexColumnToMatrix(x) 
  } else {
    rownames(x) <- rnames
    colnames(x) <- cnames 
  }
  x
}


# add names to columns and rows and do trimming
# along 1=constructs, 2=elements
#
addNamesToMatrix <- function(x, m, trim=7, along=1){
  if (!inherits(x, "repgrid")) 							    # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  if (along == 1){
    cnamesl <- getConstructNames(x)[ ,1]
    cnamesr <- getConstructNames(x)[ ,2]
    new.names <- paste(cnamesl, cnamesr, sep = " - ")
  } else {
    new.names <- getElementNames(x)
  }
  if (!is.na(trim))                               # trim constructs if prompted
    new.names <- substr(new.names, 1, trim)
  rownames(m) <- colnames(m) <- new.names         # assign new names to row and column names
  m
}

# new version using helper functions
# add names to columns and rows and do trimming
# along 1=constructs, 2=elements
#
addNamesToMatrix2 <- function(x, m, index=F, trim=7, along=1){
  if (!inherits(x, "repgrid")) 							    # check if x is repgrid object
  	stop("Object x must be of class 'repgrid'")
  if (along == 1){
    new.names <- getConstructNames2(x, index=index, trim=trim)
  } else {
    new.names <- getElementNames2(x, index=index, trim=trim)
  }
  rownames(m) <- colnames(m) <- new.names         # assign new names to row and column names
  m
}

#' add index column for neater colnames
#'
#'
#' @param x   \code{matrix} object 
#' @export
#' @keywords internal
#' @examples \dontrun{
#'    x <- matrix(1:9, 3)
#'    colnames(x) <- rep("Long names that occupiy too much space", 3)
#'    rownames(x) <- rep("Some text", 3)
#'    addIndexColumnToMatrix(x)
#' }
#'
addIndexColumnToMatrix <- function(x){ 
  if (dim(x)[1] != dim(x)[2])
    stop("works for square matrices only")
  indexes <- 1L:dim(x)[1]
  res <- cbind(indexes, x)
  colnames(res) <- c(" ", indexes)
  res
}


#' Make a histogram with steps instead of bars. Densities are used
#' for the heights.
#'
#' @title           Density histogram withs steps instead of bars
#'
#' @param vals      Numeric values to display.
#' @param breaks    Passed on to \code{hist}. 
#'                  See ?hist parameter \code{breaks} for more information.
#' @param add       Whether to add the steps to an existing plot (\code{FALSE})
#'                  or to create a new plot (default \code{add=TRUE}).  
#' @author          Mark Heckmann
#' @export
#' @keywords        internal  
#' @examples \dontrun{
#'
#'    x <- rnorm(1000) 
#'    y <- rnorm(1000, sd=.6)  
#'    stepChart(y, breaks=50)
#'    stepChart(x, add=T, breaks=50, col="red")
#' }
#'
stepChart <- function(vals, breaks="Sturges", add=FALSE, ...){
  h <- hist(vals, breaks=breaks, plot=F)
  x <- h$breaks
  y <- h$density
  x <- c(x, x[length(x)])
  y <- c(0, y, 0)
  if (add)
    points(x, y, type="s", ...) else
    plot(x, y, type="s", ...)
}


list_to_dataframe <- function(l){
  #plyr:::list_to_dataframe(l)
  do.call(rbind.data.frame, l)
}