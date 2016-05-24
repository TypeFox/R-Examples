## ----setup, include=FALSE------------------------------------------------
pkg <- 'RcppOctave'
require( pkg, character.only=TRUE )
prettyVersion <- packageVersion(pkg)
pkgtitle <- packageDescription(pkg)$Title
prettyDate <- format(Sys.Date(), "%B %e, %Y")
library(knitr)
knit_hooks$set(try = pkgmaker::hook_try)

## ----pkgmaker_preamble, echo=FALSE, results='asis'-----------------------
pkgmaker::latex_preamble()

## ----bibliofile, echo=FALSE, results='asis'------------------------------
pkgmaker::latex_bibliography(pkg)	

## ----sample_CallOctave---------------------------------------------------
.CallOctave('version')
.CallOctave('sqrt', 10)
.CallOctave('eye', 3)
.CallOctave('eye', 3, 2)

## ----ex_function, echo=FALSE, results='asis'-----------------------------
cat(readLines(system.file('scripts/ex_functions.m', package=pkg)), sep="\n")

## ----src_ex--------------------------------------------------------------
# source example function definitions from RcppOctave installation
sourceExamples('ex_functions.m')
# several functions are now defined
o_ls()

## ----call_ex-------------------------------------------------------------
# single output value
.CallOctave('fun1') 
# 3 output values 
.CallOctave('fun2')

# no output value
.CallOctave('fun_noargout', 1)
.CallOctave('fun_noargout', 'abc')

# variable number of arguments
.CallOctave('fun_varargin')
.CallOctave('fun_varargin', 1, 2, 3)

## ----argout, error = TRUE, try = TRUE------------------------------------
.CallOctave('fun_varargout')
.CallOctave('fun_varargout', argout=1)
# this should throw an error
try( .CallOctave('fun_varargout', argout=2) )

## ----argout_mod----------------------------------------------------------
# single output variable: result is S
.CallOctave('svd', matrix(1:4, 2))
# 3 output variables: results is [U,S,V]
.CallOctave('svd', matrix(1:4, 2), argout=3)
# specify output names (and therefore number of output variables)
.CallOctave('svd', matrix(1:4, 2), argout=c('U', 'S', 'V'))

## ----O_object------------------------------------------------------------
.O
.O$version()
.O$eye(3)
.O$svd(matrix(1:4,2))
# argout can still be specified
.O$svd(matrix(1:4,2), argout=3)

## ----O_object_variables, try = TRUE, error = TRUE------------------------
# define a variable
.O$myvar <- 1:5
# retrieve value
.O$myvar
# assign and retrieve new value
.O$myvar <- 10
.O$myvar
# remove 
.O$myvar <- NULL
# this should now throw an error since 'myvar' does not exist anymore
try( .O$myvar )

## ----O_object_function_call----------------------------------------------
# density of x=5 for Poisson(2)
.O$poisspdf(5, 2)
# E.g. compare with R own function
dpois(5, 2)

## ----O_object_functions--------------------------------------------------
# retrieve Octave function
f <- .O$poisspdf
f
# call (in Octave)
f(5, 2)

## ----assign, error = TRUE, try = TRUE------------------------------------
## ASSIGN
o_assign(a=1)
o_assign(a=10, b=20)
o_assign(list(a=5, b=6, aaa=7, aab=list(1,2,3)))

## GET
# get all variables
str( o_get() )
# selected variables
o_get('a')
o_get('a', 'b')
# rename on the fly
o_get(c='a', d='b')
# o_get throw an error for objects that do not exist
try( o_get('xxxxx') )
# but suggests potential matches
try( o_get('aa') )
# get a function
f <- o_get('svd')
f

## ----o_eval, error = TRUE, try = TRUE------------------------------------
# assign variable 'a'
o_eval("a=1")
o_eval("a") # or .O$a
o_eval("a=svd(rand(3))")
.O$a
# eval a list of statements
l <- o_eval("a=rand(1, 2)", "b=randn(1, 2)", "rand(1, 3)") 
l
# variables 'a' and 'b' were assigned the new values
identical(list(.O$a, .O$b), l[1:2])

# multiple statements are not supported by o_eval
try( o_eval("a=1; b=2") )
.O$a
# argument CATCH allows for recovering from errors in statement
o_eval("a=usage('ERROR: stop here')", CATCH="c=3")
.O$a
.O$c

## ----o_source------------------------------------------------------------
# clear all session 
o_clear(all=TRUE)
o_ls()

# source example file from RcppOctave installation
mfile <- system.file("scripts/ex_source.m", package='RcppOctave')
cat(readLines(mfile), sep="\n")
o_source(mfile)
# Now objects 'a', 'b', and 'c' as well as the function 'abc'
# should be defined:
o_ls(long=TRUE)
# 
o_eval("abc(2, 4, 6)")
o_eval("abc(a, b, c)")

## ----o_source_text-------------------------------------------------------
o_source(text="clear a b c; a=100; a*sin(123)")
# last statement is stored in automatic variable 'ans'
o_get('a', 'ans')

## ----o_ls----------------------------------------------------------------
o_ls()
o_ls(long=TRUE)

#clear all (variables + functions)
o_clear(all=TRUE)
o_ls()

## ----o_help, eval=FALSE--------------------------------------------------
#  o_help(std)

## ----o_doc, eval=FALSE---------------------------------------------------
#  o_doc(poisson)

## ----errors, try = TRUE--------------------------------------------------
# error
res <- try(.CallOctave('error', 'this is an error in Octave'))
geterrmessage()

# warning
res <- .CallOctave('warning', 'this is a warning in Octave')

## ----sample_svd----------------------------------------------------------

o_svd <- function(x){
	# ask for the complete decomposition
	res <- .O$svd(x, argout=c('u','d','v'))
	# reformat/reorder result
	res$d <- diag(res$d)
	res[c(2, 1, 3)]
}

# define random data
X <- matrix(runif(25), 5)

# run SVD in R
svd.R <- svd(X)
# run SVD in Octave
svd.O <- o_svd(X)
str(svd.O)
# check results
all.equal(svd.R, svd.O)
# but not exactly identical
all.equal(svd.R, svd.O, tol=10^-16)

## ----set_seed, tidy=FALSE------------------------------------------------

Rf <- function(){
	x <- matrix(runif(100), 10)
	y <- matrix(rnorm(100), 10)
	(x * y) %*% (x / y)
}

Of <- {
# define Octave function
o_source(text="
function [res] = test()
x = rand(10);
y = randn(10);
res = (x .* y) * (x ./ y);
end
")
# return the function
.O$test
}

# run both computations with a common seed  
set.seed(1234); res.R <- Rf()
set.seed(1234); res.O <- Of()
# compare results
identical(res.R, res.O)

# not seeding the second computation would give different results
set.seed(1234);
identical(Rf(), Of())


## ----news, echo=FALSE, results='asis'------------------------------------
cat(paste(readLines(system.file('NEWS', package='RcppOctave')), collapse="\n"))

## ----session_info, echo=FALSE, comment=NA--------------------------------
sessionInfo()

