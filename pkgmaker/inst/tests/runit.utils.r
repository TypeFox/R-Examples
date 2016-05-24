# Unit test for utils
# 
# Author: Renaud Gaujoux
###############################################################################

library(stringr)

test.errorCheck <- function(){
	
	f <- function(err=''){
		success <- exitCheck()
		on.exit( if(success()) cat("no error\n") else cat("with error\n") )
		
		if( err=='error' ) stop('There is an error')
		if( err=='try' ) try(stop('Catched error'), silent=TRUE)
		if( err=='tryCatch' ) tryCatch(stop('Catched error'), error = function(e){})
		
		success(1+1)
	}
	
	# without error
	out <- capture.output(res <- f())
	checkIdentical(res, 2, 'If no error: return result')
	checkIdentical(out, 'no error', 'If no error: correctly detected no error')
	
	# with error
	out <- capture.output(res <- try(f('error'), silent=TRUE))
	checkTrue( is(res, 'try-error'), 'If error: effectively throws an error')
	checkIdentical(out, 'with error', 'If error: correctly detected the error')
	
	# with try-caught error 
	out <- capture.output(res <- f('try'))
	checkIdentical( res, 2, 'If try-catched error: return result')
	checkIdentical(out, 'no error', 'If try-catched error: correctly detected no error')
	
	# with tryCatch-caught error 
	out <- capture.output(res <- f('tryCatch'))
	checkIdentical( res, 2, 'If tryCatch-catched error: return result')
	checkIdentical(out, 'no error', 'If tryCatch-catched error: correctly detected no error')
}


test.ExposeAttribute <- function(){
	
	
	x <- 1:10
	checkIdentical(ExposeAttribute(x), {attr_mode(x) <- 'rw'; x}
		, "Using ExposeAttribute() and attr_mode <- 'rw' is equivalent")
	x <- 1:10
	checkIdentical(capture.output(print(ExposeAttribute(x, a='r', b='rw'))), capture.output(print(x))
		, "Printing object with exposed attribute is identical to plain print")

	checkSet <- function(x, name, msg, ...){
		attr(x, name) <- 1
		y <- ExposeAttribute(x, ...)
		eval(parse(text=str_c('y$', name, ' <- 1')))
		attr_mode(y) <- NULL 
		checkIdentical(x, y, msg)
	}
	checkSetException <- function(x, name, msg, ...){
		y <- ExposeAttribute(x, ...)
		checkException(eval(parse(text=str_c('y$', name, ' <- 1'))), msg)
	}
	
	checkSet(x, 'a', "Set works if default")
	checkSet(x, 'a', .MODE='rw', "Set works if all args are 'rw'")
	checkSet(x, 'a', a='rw', "Set works if specified arg is 'rw'")
	checkSet(x, 'a', a='w', "Set works if specified arg is 'w'")
	checkSet(x, 'a', a='rw', b='r', "Set works if specified arg is 'rw', even if others are not")
	checkSet(x, 'ab', ab='rw', `a.*`='r', "Set works if specified arg is 'rw', even if another match is not")
	checkSetException(x, 'a', .MODE='r', "Set throws an error if access right is 'r'")
	checkSetException(x, 'a', a='r', "Set throws an error if specific access right is 'r'")
	checkSetException(x, 'a', a='', "Set throws an error if specific access right is ''")
	
	checkGet <- function(x, name, msg, ...){
		attr(x, name) <- 1
		y <- ExposeAttribute(x, ...)
		a <- eval(parse(text=str_c('y$', name)))
		checkIdentical(attr(x, name), a, msg)
	}
	checkGetException <- function(x, name, msg, ...){
		y <- ExposeAttribute(x, ...)
		checkException(eval(parse(text=str_c('y$', name))), msg)
	}
	
	checkGet(x, 'a', "Get works if default")
	checkGet(x, 'a', .MODE='rw', "Get works if all args are 'rw'")
	checkGet(x, 'a', a='rw', "Get works if specified arg is 'rw'")
	checkGet(x, 'a', a='r', "Get works if specified arg is 'r'")
	checkGet(x, 'a', a='rw', b='w', "Get works if specified arg is 'rw', even if others are not")
	checkGet(x, 'ab', ab='r', `a.*`='w', "Get works if specified arg is 'rw', even if another match is not")
	checkGetException(x, 'a', .MODE='w', "Get throws an error if access right is 'r'")
	checkGetException(x, 'a', a='w', "Get throws an error if specific access right is 'r'")
	checkGetException(x, 'a', a='', "Get throws an error if specific access right is ''")
	
	
}


test.Sys.getenv_value <- function(){
    
    on.exit( Sys.unsetenv('TOTO') )
    
    # undefined returns FALSE
    checkIdentical(Sys.getenv_value('TOTO'), FALSE, 'undefined returns FALSE')
    # raw undefined returns NA
    checkIdentical(Sys.getenv_value('TOTO', raw = TRUE), as.character(NA), 'raw undefined returns NA')
    
    Sys.setenv(TOTO='bla')
    checkIdentical(Sys.getenv_value('TOTO'), 'bla', 'defined returns value')
    
    # anything false-like returns FALSE
    Sys.setenv(TOTO='false');
    checkIdentical(Sys.getenv_value('TOTO'), FALSE, '"false" returns FALSE')
    Sys.setenv(TOTO='FALSE');
    checkIdentical(Sys.getenv_value('TOTO'), FALSE, '"FALSE" returns FALSE')
    Sys.setenv(TOTO='0');
    checkIdentical(Sys.getenv_value('TOTO'), FALSE, '"0" returns FALSE')
    
}


test.str_bs <- function(){
    
    checkIdentical(str_bs("abcd"), "abcd", "No backspace returns string unchanged")
    checkIdentical(str_bs("abcd\b"), "abc", "One backspace at the end is OK")
    checkIdentical(str_bs("\babcd"), "abcd", "One leading backspace is OK")
    checkIdentical(str_bs("abcd\b\b"), "ab", "Two backspaces at the end is OK")
    checkIdentical(str_bs("abcd\b\b\b"), "a", "Three backspaces at the end is OK")
    checkIdentical(str_bs("abcd\b\b\b\b"), "", "As many backspaces as characters at the end is OK")
    checkIdentical(str_bs("abcd\b\be"), "abe", "Backspace in the middle is OK")
}
