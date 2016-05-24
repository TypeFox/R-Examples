# Unit tests for OctaveFunctions
# 
# Author: Renaud Gaujoux
###############################################################################

library(stringr)

test.definition <- function(){
	
	# cleanup on exit
	on.exit( o_rm(all=TRUE) )
	
	.check <- function(type, ref, ...){
		msg <- function(...) str_c(type, ' [', ref, ']: ', ...)
		checkTrue(is(f <- OctaveFunction(...), 'OctaveFunction'), msg("loads without error"))
		checkIdentical(f@name, ref, msg("name of is correctly set"))
	}
	
	# built-in
	.check("Built-in function", 'randn', 'randn')
	.check("Built-in function (with arg 'name')", 'randn', fun='randn')
	# sourced
	o_source(text="function x = toto()  x=1; end")
	.check("Sourced function", 'toto', 'toto')
	.check("Sourced function (with arg 'name')", 'toto', fun='toto')
	o_rm(all=TRUE)
	# defined on the fly
	.checkruntime <- function(type, def){
		on.exit( o_rm(all=TRUE) )
		type <- str_c("Runtime function", " (", type, ")")
		.check(type, 'tata', def)
		.check(str_c(type, " (with arg 'fun')"), 'tata', fun=def)
	}
	.checkruntime("no result", "function tata()  x=1; end")
	.checkruntime("single result + no braket", "function x = tata()  x=1; end")
	.checkruntime("single result + no bracket + long name", "function xyz = tata()  x=1; end")
	.checkruntime("single result + bracket", "function [x] = tata()  x=1; end")
	.checkruntime("single result + bracket + long name", "function [xyz] = tata()  x=1; end")
	.checkruntime("multiple result", "function [x, y] = tata()  x=1; end")
	.checkruntime("multiple result + bracket + short and long names", "function [x, xyz] = tata()  x=1; end")
	
}

