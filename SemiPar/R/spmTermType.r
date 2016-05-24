########## R function: spmTermType ##########

# For determining type of term.

# Last changed: 18 JAN 2005

spmTermType <- function(term)
{
  # Determine whether or not the term is
  # a penalized function term.

  type <- NULL

  first.six.char <- substring(term,1,6)

  if( first.six.char=="offset")
     type <- "offset"

  first.two.char <- substring(term,1,2)

  if ((first.two.char!="f(")&( first.six.char!="offset"))
     type <- "lin"
  else
     {
      if ( first.six.char!="offset")
      {
        term.parts <- break.string(term,sep=",")
        if (length(term.parts)==1)     
           type <- "pen"
         else
        {
         term.parts <- break.string(term.parts[2],sep="=")
          if (length(term.parts)==1)
            type <- "krige"
         else
           type <- "pen"
        }
       }
      }

  return(type)
 }

########## End of spmTermType ##########
