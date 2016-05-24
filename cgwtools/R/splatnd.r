# function to overload "!" for one purpose only
 #this is adapted  from the sos package code for "???", credited to Duncan Murdoch.
# Basically this is a cheap sample function to show how to create a specialized unary operator that doesn't require
# parentheses for its argument.  So far as I can tell,  the only way to do this is to overload an existing function or
# operator which doesn't require parentheses.  "?" and "!" meet this requirement.
`!` <- function (e1, e2)  { 
	call <- match.call()
#  match.call breaks out each callable function in argument list (which was "??foo" for the sos package "???",
#  which allows topicExpr1 to become  a list variable w/ callable function "!"  (or "?" in sos) 

	original <- function() {
  # to call the original ? function, Duncan wrote; 
   # call[[1]] <- quote(utils::`?`)
# so I change this to: 
		call[[1]]<-quote(base::`!`)
		return(eval(call, parent.frame(2)))
	}

  # No doubt certain argument types will throw an error, and this does preclude my ever having an actual
  # variable called "newdev" (or at least trying to create the actual NOT of it) 
  # Interesting: when using "!" for its real purpose, e.g. " !(bar2 %in% bar1) ", converting to
  # character will create multiple elements, so the "collapse" is critical here
	switch(paste(as.character(call[[2]]),sep='',collapse='') ,
	'newdev' = dev.new(width=4.5, height= 4.5, restoreConsole=T),
	'qapla' = cat('batlh tIn chav\n'),
	return(original()) )
	
}

