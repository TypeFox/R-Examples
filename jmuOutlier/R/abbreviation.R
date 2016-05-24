abbreviation <-
function(x, choices) {
  # Returns a value in 'choices' specified by 'x', which may be an abbreviation.
  # If no such abbreviation exists, then the original value of 'x' is returned.
  # 'x': A character string, and consists of some or all letters in a value in \code{choices} or may equal \code{choices}.
  # 'choices': A vector of character strings.
  # Note:  \code{x} is typically a scalar.  However, if \code{x} is a vector, then the first value is used.
  # Examples:
  #     x1 = "two";  x2 = "l";   x3 = "gr";   x4 = "greater";  x5 = "NotInChoices"
  #     choices = c("two.sided", "less", "greater")
  #     abbreviation( x1, choices ) ; abbreviation( x2, choices ) ; abbreviation( x3, choices ) 
  #     abbreviation( x4, choices ) ; abbreviation( x5, choices )
  if ( !is.character(x) )   stop("'x' must be a character string.")
  if ( !is.character(choices) )   stop("'choices' must be a character string.")
  answer <- x[1]
  for (i in 1:length(choices))  { if (!is.na(charmatch(x[1], choices[i])))  answer <- choices[i] }
  return(answer)
}
