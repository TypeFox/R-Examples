#{{{ A test of inlinedocs with some difficult items
'@
To compile the package, make:
<code R>
library(inlinedocs)
package.skeleton.dx("~/Documents/Pgm/Rforge/inlinedocs/pkg/inlinedocs/scratch/inlinedocsTest", namespace = TRUE)
</code>
<<>>='
#}}}


"%!in%" <- function ( # Look which elements are absent from a table
### Check which ones of the tested values are absent from a reference table.
### This is the opposite as %in% operator in base package.
x,		##<< a vector of values to check
table	##<< a vector of items of the same type as 'x', to be checked
){
	return(!x %in% table)
	### A vector of logical of the same length as 'x', indicating which items
	### are not in 'table'
}

"[<-.ex" <- function ( # Replace an item in an 'ex' object
x,			##<< an 'ex' object
which = 1,	##<< which element of the 'ex' object to replace
value		##<< the value to use for replacement
){
	# This is not much useful, but it is used to test a special case of usage!
	x <- as.character(x)
	x[which] <- value
	class(x) <- c("ex", "character")
	return(x)
	### 'x' is returned, with 'value' replacing the 'which'th element
}

# Note: to embed examples in objects as 'ex' attribute:
"ex<-" <- function (  # Embed examples code in an object
x,		##<< an object
value	##<< a character string containing code for examples
){
	# Add 'value' as a 'ex' attribute to 'x', after coercing to an 'ex' object
	attr(x, "ex") <- structure(as.character(value), class = c("ex", "character"))
	return(x)
	### The initial object 'x' with attribute 'ex' added or changed to 'value'
}

"is.ex" <- function (	# Test if an object is of 'ex' class, or possesses an 'ex' attribute
x		##<< an object to test
){
	return(inherits(x, 'ex') || !is.null(attr(x, "ex")))
	### Return \code{TRUE} if 'x' inherits of class 'ex' or has an 'ex' attribute,
	### \code{FALSE} otherwise.
}

"ex" <- function (	# Extract and run the 'ex' attribute from an object (or run the 'ex' object itself)
x,				##<< an object possibly containing an 'ex' attribute, or of class 'ex'
runIt = TRUE,	##<< do we run the examples?
...				##<< further arguments passed to source()
){
	if (inherits(x, 'ex')) {
		res <- x
	} else res <- attr(x, "ex")
	if (is.null(res)) {
		cat("No examples found\n")
		return(invisible(FALSE))
	} else {
		res <- structure(as.character(res), class = c("ex", "character"))
		if (isTRUE(runIt)) {
			txt <- textConnection(res)
			on.exit(close(txt))
			cat("Examples for '", deparse(substitute(x)), "'\n", sep = "")
			source(txt, echo = TRUE, ...)
		}
		return(invisible(res))
	}
	### The extracted 'ex' object
}

"print.ex" <- function (	# Print method of an 'ex' object
x,		##<< an object that has a non empty 'ex' attribute
...		##<< further arguments to pass to print()
){
	if (!is.ex(x)) stop("'x' must be an 'ex' object")
	cat("Examples:\n")
	cat(x, "\n")
	return(invisible(x))
}

"Combine" <- function (x, y, ...) # Combine lists or character strings
	UseMethod("Combine")

"Combine.character" <- function ( # Combine character strings by pasting them together
x,	##<< a first character vector to combine
y,	##<< a second character vector to combine
...	##<< further arguments (not used yet)
){
	res <- paste(x, y, sep = "\n")
	return(res)
	### The comment of the return value
	### on two lines...
	#{{{examples
	a <- c("A", "B", "C")
	b <- c("A", "C", "D")
	Combine(a, b)
	#}}}examples
}

"Combine.list" <- structure(function( # Combine lists by adding elements or adding to existing elements
x,		##<< a first list to combine
y,		##<< a second list to combine
		##   with the previous one 
...		##<< further arguments (not used yet)
#{{{doc
#S3method<<Combine
##description<< This is a more extensive description of what the function
##				is supposed to do... It typically spans on several lines
##				and the intelligent formatting of Komodo is working very
##				nicely in this case (use Shift+Return)
##
##note<< This is just to test inlinedocs extended features inside SciViews
##		 Komodo, including formattings like \eqn{n}th power.
##
##references<< \url{http://inlinedocs.r-forge.r-project.org}
##
##seealso<< \code{\link{combine}}
#}}}doc
){
	toadd <- !names(y) %in% names(x)
	toup <- names(y)[names(y) %in% names(x)]
	x[names(y)[toadd]] <- y[toadd]
	for (up in toup)
		x[[up]] <- combine(x[[up]], y[[up]])
	return(x)
	### A list that combines both 'x' and 'y' lists
},
# Here is how to attach examples code to an object
ex = function () {
	# Simple example of combining two lists
	a <- list(x = "hello", y = "world")
	b <- list(x = "John", z = "everybody")
	Combine(a, b)
})
