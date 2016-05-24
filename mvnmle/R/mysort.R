mysort<-function(x)
  {
    # Sorts rows and cols of incoming dataframe x into/
    # an order for which it is easier to write the likelihood function
    nvars<-ncol(x)
    powers<-as.integer(2^((nvars-1):0))
    binrep<-ifelse(is.na(x),0,1)
    decrep<-binrep %*% powers
    sorted<-x[order(decrep),]
    decrep<-decrep[order(decrep)]

    list(sorted.data=sorted, freq=as.vector(table(decrep)))
  }
