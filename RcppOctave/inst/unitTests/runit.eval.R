# Unit tests for eval functions: assign, get, etc...
# 
# Author: Renaud Gaujoux
# Creation: 17 Nov 2011
###############################################################################


#' Unit test for o_eval
test.o_eval <- function(){
	
	o_clear(all=TRUE)
	val <- 10
	checkIdentical(o_eval('a=10'), val, "returns result")
	checkIdentical(o_eval('a'), val, "value is effectively assigned and retrieved")
	
	# unlist
	checkIdentical(o_eval('a', unlist=FALSE), list(val), "variables is retrieved as unamed list if unlist=FALSE")
	
	# errors
	checkException(o_eval('a=3; b=5;'), "Error if multiple statements")
	checkException(o_eval('b'), "Error with no CATCH argument: throws an error")
	
	# CATCH
	checkIdentical(o_eval('b', CATCH='c=2*10'), 2*val, "No error with no CATCH argument and value of CATCH is returned")
	checkIdentical(o_eval('c'), 2*val, "Value of CATCH was evaluated")
	
	# list of statement
	checkIdentical(o_eval('a=1', 'b=5', 'c=0.8'), list(1,5,0.8), "Multiple arguments are evaluated as multiple statements")
	val <- list(10,50,0.08)
	checkIdentical(o_eval(list('a=10', 'b=50', 'c=0.08')), val, "Single list argument is evaluated as multiple statements: returns values")
	checkIdentical(o_eval(list('a', 'b', 'c')), val, "Single list argument is evaluated as multiple statements: variable assignment are effective")
	nam <- c('A','B', 'C')
	checkIdentical(o_eval(setNames(list('a', 'b', 'c'), nam)), setNames(val, nam), "Single named list argument: evaluation + names")
	
}

#' Unit test for o_get
test.o_get <- function(){
	
	o_clear(all=TRUE)
	o_eval('who')
	
	checkIdentical(o_get(), list(), "No argument + empty context: returns empty list")
	checkIdentical(o_get(rm.ans=FALSE), list(ans=o_eval('ans')), "No argument + empty context + rm.ans=FALSE: returns 'ans' value in a list")
	checkIdentical(o_get('ans'), o_eval('ans'), "o_get('ans'): returns 'ans' value")
	checkIdentical(o_get('ans', rm.ans=TRUE), list(), "o_get('ans', rm.ans=TRUE): returns empty list")
	
	o_clear()
	o_eval('who')
	l <- list(a=1)
	o_load(l)
	checkIdentical(o_get(), l, "No argument + single variable in context: returns single variable in list")
	checkIdentical(o_get(unlist=TRUE), l$a, "No argument + single variable in context + unlist: returns single variable value")
	checkIdentical(o_get('a'), l$a, "Single quoted argument + single variable in context: returns single variable value")
	checkIdentical(o_get('a', unlist=FALSE), l, "Single unquoted argument + single variable in context + unlist: returns single variable in list")
	
	o_clear()
	o_eval('who')
	l <- list(b=1, c=3, d=matrix(1:9, 3))
	o_load(l)
	
	checkIdentical(o_get(), l, "No argument + multiple variable in context: returns all variables in list")
	checkIdentical(o_get(rm.ans=FALSE), c(list(ans=o_get('ans')), l), "No argument + rm.ans=FALSE: returns all variables in list including 'ans'")
	checkIdentical(o_get('b'), l$b, "Single quoted argument: returns single variable value")
	checkIdentical(o_get('b', unlist=FALSE), l['b'], "Single unquoted argument: returns single variable value in list")
		
	checkIdentical(o_get(X='b'), list(X=l$b), "Single named argument: returns renamed variable in list")
	checkIdentical(o_get(X='b', unlist=TRUE), l$b, "Single named argument: returns variable value")
	checkIdentical(o_get(X='b', Y='c'), list(X=l$b, Y=l$c), "Two named arguments: returns renamed variables")
	checkIdentical(o_get(X='b', 'c', 'd'), list(X=l$b, c=l$c, d=l$d), "Mix of named/unnamed arguments: returns renamed/named variable list")
	
	# check using $ and [[
	checkIdentical(o_get('b'), .O$b, 'Binding of .O$ to o_get works')
	checkIdentical(o_get('b'), .O[['b']], 'Binding of [[ to o_get works')
	
	f <- o_get('svd')
	checkTrue( is(f, 'OctaveFunction'), "Get of a function returns an OctaveFunction object")
	checkIdentical( f@name, 'svd', "Get of a function returns the correct OctaveFunction object")
	
	f <- o_get(a='svd')
	checkTrue( is(f$a, 'OctaveFunction'), "Get of a function with named argument returns a list with the OctaveFunction object in named element")
	checkIdentical( f$a@name, 'svd', "Get of a function with named argument returns the correct OctaveFunction object in list")
	
	# Errors in matching
	o_clear(all=TRUE)
	o_assign(a=1, aaaa=2, aab=3)
	checkIdentical(o_get('a'), 1, "Ok if multiple partial matches but one exact match")
#	checkIdentical(o_get('a', exact=TRUE), 1, "Ok if multiple partial matches but one exact match (arg exatc=TRUE)")
	checkException(o_get('aa'), "Error if multiple matches but no exact match")
#	checkException(o_get('aa', exact=TRUE), "Error if multiple matches but no exact match (exact = TRUE)")
#	checkIdentical(o_get('aaa'), 2, "Ok if unique partial match but no exact match")
#	checkException(o_get('aaa', exact=TRUE), "Error if unique partial match but no exact match and exact=TRUE ")
	checkException(o_get('aaa'), "Error if unique partial match but no exact match")
	
}


test.o_source <- function(){
	
	o_clear(all=TRUE)
	# text argument
	ref <- list(a=1, b=3, c=as.numeric(1:3))
	txt <- c("a=1", "b=3", "c=[1 2 3]")
	o_source(text=txt)
	checkIdentical(o_get('a','b','c'), ref, "'text' argument works for evaluating multiple statements passed as a character vector")
	o_clear()
	o_source(text=paste(txt, collapse=";"))
	checkIdentical(o_get('a','b','c'), ref, "'text' argument works for evaluating multiple statements passed as a single character string")
}

test.assign <- function(){

	o_clear(all=TRUE)
	
	## assign a variable as named arguments
	o_assign('a', 50)
	checkIdentical(o_get('a'), 50, "Assign with two unamed arguments: character name + value works")
	
	## assign a variable as named arguments
	val <- list(a=1, b=2, c=matrix(1:9, 3))
	do.call(o_assign, val)
	checkIdentical(o_get(), val, "Assign with named argument works")
	
	## assign a variable for each element in a list
	x <- list(a=10, b=20, c=matrix(1:9 * 10, 3))
	o_assign(x)
	checkIdentical(o_get(), x, "Assign with single named list argument works")
	
	## assign the content of an environment
	x <- list(a=100, b=200, c=matrix(1:9 * 100, 3))
	e <- list2env(setNames(x, names(x)))
	o_assign(e)
	checkIdentical(o_get(), x, "Assign with from an environment")
}

test.redirection <- function(){
    
    if( .Platform$OS.type == 'windows' ) DEACTIVATED("Redirection does not currently work on Windows")
    
    # Output
    out <- 'This is some Octave text output'
    checkIdentical(capture.output(dummy <- .CallOctave('printf', out)), out, "Octave output is correctly captured")
    checkIdentical(capture.output(dummy <- .CallOctave('printf', out, buffer.std = 1)), out, "Octave output is correctly captured (buffer.std = 1)")
    # no buffering
    checkTrue(identical(capture.output(dummy <- .CallOctave('printf', out, buffer.std = -1)), out), "Octave output is still captured if stdout is not buffered")
    
    # Errors
    oerr <- 'This is an Octave error indeed!!!'
    checkException(.CallOctave('error', oerr), "R error is raised by Octave error")
    checkTrue(grepl(oerr, geterrmessage(), fixed = TRUE), "Octave error message is correctly passed to R")
    checkException(.CallOctave('error', oerr), "R error is raised by Octave error (buffer.std = 2)")
    checkTrue(grepl(oerr, geterrmessage(), fixed = TRUE), "Octave error message is correctly passed to R (buffer.std = 2)")
    # no buffering
    checkException(.CallOctave('error', oerr, buffer.std = -2), "R error is raised by Octave error even when stderr is not buffered")
    checkTrue(grepl(oerr, geterrmessage(), fixed = TRUE), "Octave error message is still passed to R if stderr is not buffered")
    
    # Warnings
    checkWarning(.CallOctave('warning', 'aaaa'), "aaaa", "Default call buffer warnings")
    checkWarning(.CallOctave('warning', 'aaaa', buffer.std = 2), "aaaa", "Default call buffer warnings (buffer.std = 2)")
    checkWarning(.CallOctave('warning', 'aaaa', buffer.std=-2)
                                               , "aaaa", "Warning are still detected if not buffering stderr")
    checkIdentical(capture.output(dummy <- .CallOctave('warning', 'aaaa')), character(), "Octave warning printed directly (not in stdout) if buffering is disabled (buffer.std=0)")
    
}
