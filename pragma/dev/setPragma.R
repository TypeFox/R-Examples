# -------------------------------------------------------------------------
# setPragma:
#   Creates a keyword that can then be used to implement simple 
#   functionality with a clean syntax
#
#   R is (mostly) a functional language.  Any action taken is made through
#   function calls. This requires all of the calls to use parens even if 
#   no arguments are necessary.  Very often a function has no arguments
#   and are used to set some functionality of the program.  These are 
#   sometimes called pragmas or directives and specify some behaviours of
#   the program -- usually how to compile the program.  Of course, R is
#   not compiled.  It is completely non-essential, but having the ability
#   to have PRAGMAS is a nice feature and makes the code more readable.
#    
#    > AUTOHELP     # instead of autohelp()
# 
#   There should be some guidance on syntaz such as keywords 
#   should be in ALL CAPS.  
#
#   setPragma( function )
#    - Defines an anonymous class
#    - Defines a keyword an object of this class
#    - Defines a show method.
#
#   TODO:
#    - sealed=TRUE : seal the PRAGMA class?
#
# -------------------------------------------------------------------------

# library(formula.tools) 

setPragma <- function( fun ) {

  name <- as.character( lhs( last.call()  ) )
  name <- paste( "PRAGMA", name, sep="." ) 
  # TEST IF ALREADY PRESENT

  # SET CLASS
  setClass( name, "NULL", where=globalenv() ) 


  # CONSTRUCT FUNCTION FOR SHOW METHOD
  x <- function(object) {} 
  body(x) <- body(fun) 

  # DEFINE THE ACTION FOR THE KEYWORD 
  setMethod( "show", name, x, where=globalenv() )  

  # return(x) 
  # RETURN THE KEYWORD 
  return( new( name ) ) 

}
 

# FOO <- setPragma( function() cat( "FOO WORKS\n" )  ) 
# export(xdx)
# FOO 





