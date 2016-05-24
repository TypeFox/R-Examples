"labcon" <- function(x)
  {
      ## Verifications
      if (!inherits(x, "asc"))
          stop("should be an object of class asc")
      y<-x

      ## rajfond adds empty lines and columns on the borders of a map
      rajfond <- function(x) {
          nr <- nrow(x)
          nc <- ncol(x)
          f <- rep(0, nr)
          x <- cbind(f, x, f)
          f <- rep(0, nc + 2)
          x <- rbind(f, x, f)
      }

      ## The map is transformed so that it takes either
      ## the value 0 (NA) or 1 (mapped value)
      x[!is.na(x)] <- 1
      x[is.na(x)] <- 0
      x <- rajfond(x)

      ## sequential labelling of connex components
      ## with the C function "seqeticorr"
      toto <- .C("seqeticorr", as.double(t(x)), as.integer(nrow(x)),
                 as.integer(ncol(x)), PACKAGE="adehabitat")

      ## output
      etiquete <- matrix(toto[[1]], nrow = nrow(x), byrow = TRUE)
      ## and we delete the empty lines and columns added
      etiquete <- etiquete[-c(1, nrow(etiquete)), -c(1, ncol(etiquete))]
      etiquete[etiquete==0]<-NA
      s<-getascattr(y, etiquete)
      attr(s, "type")<-"factor"
      attr(s, "levels")<-as.character(1:nlevels(factor(etiquete)))
      return(s)
  }

