"lowres.kasc" <- function(x, np=2, ...)
  {
    if (!inherits(x, "kasc"))
      stop("x sould be of class \"kasc\"")
    so <- list()
    for (i in names(x)) {
      so[[i]]<-lowres.asc(getkasc(x, i), np)
    }
    x<-as.kasc(so)
    return(x)
  }

