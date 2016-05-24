
strmacro <- function(..., expr, strexpr)
{
  if(!missing(expr))
    strexpr <- deparse(substitute(expr))
  
  a <- substitute(list(...))[-1]

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
      else
        {
          a[[i]] <- a[[i]]
        }
    }
  names(a) <- nn
  a <- as.list(a)

  ## this is where the work is done
  ff <- 
    function(...)
      {
        ## build replacement list
        reptab <- a # copy defaults first
        reptab$"..." <- NULL
        
        args <- match.call(expand.dots=TRUE)[-1]
                          
        for(item in names(args))
          reptab[[item]] <- args[[item]]
        
        ## do the replacements
        body <- strexpr
        for(i in 1:length(reptab))
          {
            pattern <- paste("\\b",
                             names(reptab)[i],
                             "\\b",sep='')
            
            value <- reptab[[i]]
            if(missing(value))
              value <- ""
            
            body <- gsub(pattern,
                         value,
                         body)
          }

        fun <- parse(text=body)
        eval(fun, parent.frame())

        
      }
  
  
  
  ## add the argument list
  formals(ff) <- a
  
  ## create a fake source attribute
  mm <- match.call()
  mm$expr <- NULL
  mm[[1]] <- as.name("macro")
  attr(ff, "source") <- c(deparse(mm), strexpr)
  
  ## return the 'macro'
  ff
}




