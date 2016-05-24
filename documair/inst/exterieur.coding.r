
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sum4c <- function (a,b)
#TITLE example of calling  C and Fortran functions
#DESCRIPTION 
# An example of an \samp{R} function calling a \samp{C} function.
#DETAILS
# The documentation must be given within the \samp{R} function,
# even there is no need to explain that a \samp{C} function was used.
#ALIAS sum4 sum4c
#KEYWORDS 
#INPUTS
#{a}<< First integer>>
#{b}<< Second integer>>
#[INPUTS]
#VALUE
# The sum of the two integers.
#EXAMPLE
# sum4c(1,1)
# sum4c(sum4c(1,1),sum4c(2,2))
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_05_23
#REVISED 14_05_23
#--------------------------------------------
{
  # just calling the C function
  toto <- 3;
  res <- .C("sumc",as.integer(a),as.integer(b),retour=as.double(toto));
  # returning
  res$retour;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
sum4f <- function (a,b)
#TITLE example of calling a Fortran function
#DESCRIPTION 
# An example of an \samp{R} function calling a \samp{Fortran} function.
#DETAILS
# The documentation must be given within the \samp{R} function,
# even there is no need to explain that a \samp{Fortran} function was used.
#ALIAS
# sum4f
#KEYWORDS 
#INPUTS
#{a}<< First integer>>
#{b}<< Second integer>>
#[INPUTS]
#VALUE
# The sum of the two integers. 
#EXAMPLE
# sum4f(1,1)
# sum4f(sum4c(1,1),sum4c(2,2))
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 14_05_23
#REVISED 14_05_23
#--------------------------------------------
{
  # just calling the Fortran function
  uu <- 67;
  res <- .Fortran("sumf",as.integer(a),as.integer(b),voir=as.double(uu));
  # returning
  res$voir;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
