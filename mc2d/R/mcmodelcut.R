#<<BEGIN>>
mcmodelcut <- function(x, is.expr=FALSE)
#ISALIAS evalmccut
#--------------------------------------------
#
{

  if(!is.expr) x <- as.expression(substitute(x))
  if(!is.expression(x)) stop("x can not be evaluate as an expression")

  nbexpr <- length(x[[1]])
  if(nbexpr != 4) stop("The expression should include three blocks.")

  cat("The following expression will be evaluated only once :\n")
  cat(deparse(x[[1]][[2]]),sep="\n")

  last <- x[[1]][[3]][[length(x[[1]][[3]])]][[3]]
  lastcall1 <- substr(deparse(last,width.cutoff = 500), 1, 3)
  if (lastcall1 != "mc(") warning("The last call should be 'mymc <- mc(...)'")

  nom <- as.character(x[[1]][[3]][[length(x[[1]][[3]])]][[2]])                  # name of the mc object in the third block
  cat("The mc object is named: ",nom,"\n")

  class(x) <- "mcmodelcut"
  return(invisible(x))
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
