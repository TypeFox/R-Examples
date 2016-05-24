ci.dose.radir <-
  function(object, cr=0.95)
  {
    if (class(object)!="dose.radir") stop("Wrong object")
    if (cr<=0 | cr>=1) stop("Wrong confidence region")
    
    cdaux <- approxfun(object[[2]], object[[1]])
    cd <- function(x){
      ifelse(is.na(cdaux(x)), 0, cdaux(x))}
    imod <- which.max(object[[1]])
    y <- seq(.001, cd(object[[2]][imod]) - .001, .001)
    lowd <- ifelse(0<object[[2]][imod],uniroot(function(x) y[1] - cd(x), c(0, object[[2]][imod]), extendInt="yes")$root,0)
    uppd <- ifelse(class(try(uniroot(function(x) y[1] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes"), silent=T))!="try-error",
                   uniroot(function(x) y[1] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes")$root, max(object[[2]]))
    i <- 1
    while(integrate(Vectorize(function(x) cd(x)), lowd, uppd)$value > cr)
    {
      i <- i + 1
      lowd <- ifelse(0<object[[2]][imod],uniroot(function(x) y[i] - cd(x), c(0, object[[2]][imod]), extendInt="yes")$root,0)
      uppd <- ifelse(class(try(uniroot(function(x) y[i] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes"), silent=T))!="try-error",
                     uniroot(function(x) y[i] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes")$root, max(object[[2]])) 
    }
    lowd <- ifelse(0<object[[2]][imod],uniroot(function(x) y[i-1] - cd(x), c(0, object[[2]][imod]), extendInt="yes")$root,0)
    if (lowd < 0) lowd <- 0
    uppd <- ifelse(class(try(uniroot(function(x) y[i-1] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes"), silent=T))!="try-error",
                   uniroot(function(x) y[i-1] - cd(x), c(object[[2]][imod], max(object[[2]])), extendInt="yes")$root, max(object[[2]]))
    return(c(lowd, uppd))
  }