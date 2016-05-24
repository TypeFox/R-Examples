"gvlma" <-
function(x, data, alphalevel = 0.05, timeseq, ...)
  {
    mc <- match.call()
    extras <- match.call(expand.dots = FALSE)$...
    
    if (class(x) == "formula")
      {
        if (missing(timeseq)) timeseq = 1:nrow(data)
        calllist <- call(name = "gvlma.form",
                         formula = x,
                         data = as.name(deparse(substitute(data))),
                         alphalevel = alphalevel,
                         timeseq = timeseq)
        if (length(extras)>0)
          {
            calllist <- c(as.list(calllist), extras)
            calllist <- as.call(calllist)
          }
        z <- eval(calllist, parent.frame())
      }
    else
      z <- gvlma.lm(x, alphalevel, timeseq)
    z$GlobalTest$call <- mc
    z
  }

