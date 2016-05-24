"kasc2df" <- function(x, var=names(x))
  {
      ## Verifications
      if (!inherits(x, "kasc")) stop("Non convenient data type")

      ## Bases
      w<-data.frame(x[var])
      index<-c(1:nrow(w))

      ## abenner returns TRUE if a vector contains no missing values
      abenner<-function(x){
          if (any(is.na(x))) {
              return(FALSE)
          } else {
              return(TRUE)
          }
      }

      ## deletes all the rows containing NAs
      cons<-apply(w, 1, abenner)
      indcons<-index[cons]
      wcons<-data.frame(w[cons,])

      ## and return the index
      output<-list(index=indcons, tab=wcons)
  }

