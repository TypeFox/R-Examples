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

.result <- 
 list(silly.example = list(definition = "silly.example <- function\n### this function does nothing in particular and does it very well\n(\n ##title<<Simple function arguments\n first,         ##<< the first argument with a multi-line description\n ## which I really have to put here rather than explaining in the details.\n second=        ##<< the second argument with a list default value\n ## and descriptions of each of the elements\n ##describe<<\n list(this=\"that\", ##<< whichness\n      the=\"other\", ##<< of the\n      rhubarb=\"stew\", ##<< why\n      foo=\"bar\"),\n ##end<<\n third ##<< an argument that does nothing\n )\n{\n  ##description<<why should I add to description?\n  ##details<<\n  ## if second is TRUE then first is returned\n  if ( second ){\n    ##alias<<Long silly alias\n    res <- first\n  } else {\n    ##details<<\n    ## if second is not TRUE then a list is returned\n    ##describe<<The contents of the list are:\n    res <- list(x=7, ##<< x coordinate\n                z= ##<< z describes everything else\n                ##describe<<\n                list(colour=green, ##<< colour of line\n                     width=2),     ##<< width of line\n                ##end<<\n                ## and this line should get into documentation for z\n                y=10)##<< y coordinate\n  }\n  ##note<< a note\n  ##references<< a reference\n  ##seealso<< \\code{\\link{Silly-class}}\n  ##keyword<<documentation utilities\n  invisible(res)\n### invisible something not unrelated to first\n}",
        format="",
     description = "this function does nothing in particular and does it very well\nwhy should I add to description?",  
     value = "invisible something not unrelated to first", title = "Simple function arguments",  
     `item{first}` = "the first argument with a multi-line description\nwhich I really have to put here rather than explaining in the details.",  
     `item{second}` = "the second argument with a list default value\nand descriptions of each of the elements\\describe{\n\\item{this}{whichness}\n\\item{the}{of the}\n\\item{rhubarb}{why}\n}",  
     `item{third}` = "an argument that does nothing", details = "if second is TRUE then first is returned\n\nif second is not TRUE then a list is returned\n\nThe contents of the list are:\\describe{\n\\item{x}{x coordinate}\n\\item{z}{z describes everything else\\describe{\n\\item{colour}{colour of line}\n\\item{width}{width of line}\n}\nand this line should get into documentation for z}\n\\item{y}{y coordinate}\n}",  
     alias = "Long silly alias", note = "a note", references = "a reference",  
     seealso = "\\code{\\link{Silly-class}}", keyword = "documentation}\n\\keyword{utilities"),  
     `Silly-class` = list(title = "S4 classes can be documented as well",
       format="",
         description = "The Silly class does nothing much either",  
         details = "Put what you like in documentation details,\nbut ideally reference construction methods.",  
         `item{forwards}` = "forward operation", `item{reverse}` = "how to go backward",  
         `item{crashes}` = "how many crashes", `section{Objects from the Class}` = "Put what you like in documentation details,\nbut ideally reference construction methods.",  
         seealso = "", alias = "Silly")) 
