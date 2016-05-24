".managNAkasc" <- function(x)
  {
      ## Verifications
      if (!inherits(x,"kasc")) stop("non convenient data")
      class(x)<-"data.frame"

      ## We only keep the mapped pixels for ALL variables
      tmpy<-is.na(x)
      tmp<-apply(tmpy, 1, function(x) sum(as.numeric(x)))
      x[tmp!=0,]<-rep(NA, ncol(x))
      class(x)<-c("kasc", "data.frame")
      return(x)
  }

