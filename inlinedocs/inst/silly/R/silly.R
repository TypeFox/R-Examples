silly.example <- function
### this function does nothing in particular and does it very well
(
 ##title<<Simple function arguments
 first,         ##<< the first argument with a multi-line description
 ## which I really have to put here rather than explaining in the details.
 second=        ##<< the second argument with a list default value
 ## and descriptions of each of the elements
 ##describe<<
 list(this="that", ##<< whichness
      the="other", ##<< of the
      rhubarb="stew", ##<< why
      foo="bar"),
 ##end<<
 third ##<< an argument that does nothing
 )
{
  ##description<<why should I add to description?
  ##details<<
  ## if second is TRUE then first is returned
  if ( second ){
    ##alias<<Long silly alias
    res <- first
  } else {
    ##details<<
    ## if second is not TRUE then a list is returned
    ##describe<<The contents of the list are:
    res <- list(x=7, ##<< x coordinate
                z= ##<< z describes everything else
                ##describe<<
                list(colour=green, ##<< colour of line
                     width=2),     ##<< width of line
                ##end<<
                ## and this line should get into documentation for z
                y=10)##<< y coordinate
  }
  ##note<< a note
  ##references<< a reference
  ##seealso<< \code{\link{Silly-class}}
  ##keyword<<documentation utilities
  invisible(res)
### invisible something not unrelated to first
}

setClass("Silly", # S4 classes can be documented as well
### The Silly class does nothing much either
         ##details<< Put what you like in documentation details,
         ## but ideally reference construction methods.
         representation(forwards="function", ##<< forward operation
                        reverse="function", ##<< how to go backward
                        crashes="integer") ##<< how many crashes
         ) ##<< this comment is ignored as it is outside setClass expression

## creates "show" generic function. Documentation of this not yet supported.
##setMethod("show","Silly",function(object){cat("crashed ",object@crashes," times\n")})

