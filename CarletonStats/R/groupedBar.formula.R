groupedBar.formula <-
function(formula, data = parent.frame(), subset,  ...)
{
   
   if (length(formula) > 3L)
        stop("'formula' incorrect")
  
    m <- match.call(expand.dots = FALSE)
    if (is.matrix(eval(m$data, parent.frame()))) 
        m$data <- as.data.frame(data)
    m[[1L]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    varnames <- names(mf)
    
    nmiss <- length(attr(mf, "na.action"))
    if (nmiss > 0)
     cat("\n ", nmiss, "observation(s) removed due to missing values.\n")
     
     if (length(formula) == 2L) {
        do.call("groupedBar", c(list(x = unlist(mf)), list(...)))
  
      } else {
  
     response <- attr(attr(mf, "terms"), "response")
     y <- mf[[response]]
     g <- mf[[-response]]    
     names(g) <- varnames[2]
       
     y <- do.call("groupedBar", c(list(y, g), list(...)))
    }#end else
}
