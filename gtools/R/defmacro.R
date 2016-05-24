## Code from
##
## @Article{Rnews:Lumley:2001,
##  author       = {Thomas Lumley},
##  title	       = {Programmer's Niche: Macros in {R}},
##  journal      = {R News},
##  year	       = 2001,
##  volume       = 1,
##  number       = 3,
##  pages	       = {11--13},
##  month	       = {September},
##  url	       = {http://CRAN.R-project.org/doc/Rnews/}
##}
defmacro <- function(..., expr) #, DOTS=FALSE)
{
  expr <- substitute(expr)
  a <- substitute(list(...))[-1]

  ## process the argument list
  nn <- names(a)
  if (is.null(nn))
    nn <- rep("", length(a))
  for(i in 1:length(a))
    {
      if (nn[i] == "")
        {
          nn[i] <- paste(a[[i]])
          msg <- paste(a[[i]], "not supplied")
          a[[i]] <- substitute(stop(foo),
                               list(foo = msg))
        }
      if (nn[i] == "DOTS")
        {
          nn[i] <- "..."
          a[[i]] <- formals(function(...){})[[1]]
        }
    }
  names(a) <- nn
  a <- as.list(a)

  ## this is where the work is done
  ff <- eval(substitute(
                        function()
                        {
                          tmp <- substitute(body)
                          eval(tmp, parent.frame())
                        },
                        list(body = expr)))
  
  ## add the argument list
  formals(ff) <- a
  
  ## create a fake source attribute
  mm <- match.call()
  mm$expr <- NULL
  mm[[1]] <- as.name("macro")
  attr(ff, "source") <- c(deparse(mm),
                          deparse(expr))
  
  ## return the 'macro'
  ff
}
