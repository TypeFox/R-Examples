"colasc" <- function(x, ...)
  {
      ## Verifications
      if (!inherits(x, "asc"))
          stop("Should be an \"asc\" object")

      ## the list of colors with levels of maps
      l<-list(...)
      n<-names(l)
      i<-NA

      ## levels of the map
      tc<-levels(x)

      ## Verifications that the names of "..." correspond to the
      ## levels of x
      for (i in n) {
          if (!any(tc==i))
              stop(paste(i, "is not a valid level of the factors"))
      }

      ## creates the vector of colors
      coul<-0
      for (i in 1:length(tc)) {
          u<-tc[i]
          coul[i]<-l[[u]]
      }
      return(coul)
  }

