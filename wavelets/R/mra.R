mra <- function(X, filter="la8", n.levels, boundary="periodic", fast=TRUE, method="dwt"){

  # error checking
  if(is.na(match(method, c("dwt","modwt"))))
    stop("Invalid argument value for 'method'")

  # initializations
  J <- n.levels
  details <- as.list(rep(NA, length=J))
  names(details) <- lapply(1:J, function(j, x){names(x)[j] <- paste("D",j,sep="")}, x=details)
  smooths <- as.list(rep(NA, length=J))
  names(smooths) <- lapply(1:J, function(j, x){names(x)[j] <- paste("S",j,sep="")}, x=smooths)

  # compute wavelet transform and other necessary values
  if(method == "dwt"){
    wt <- dwt(X, filter, n.levels, boundary, fast)
  } else wt <- modwt(X, filter, n.levels, boundary, fast)
  filter <- wt@filter
  n.series <- dim(wt@W[[1]])[2]
  N <- dim(wt@series)[1]

  # compute the details and smooths
  for(j in 1:J){
    Wj <- wt@W[[j]]
    Vj <- wt@V[[j]]
    DWj <- Wj
    SWj <- matrix(rep(0, length=length(Wj)), ncol=n.series)
    DVj <- matrix(rep(0, length=length(Wj)), ncol=n.series)
    SVj <- Vj

    # compute the details and smooths for level j
    for(k in j:1){
      if(fast){
        if(method == "dwt"){
          Vout <- rep(0, length=2*dim(DWj)[1])
          DVj <- sapply(1:n.series,
                        function(i,w,v,f,Vout){
                          return(.C("dwt_backward", as.double(w[,i]),
                                    as.double(v[,i]), as.integer(length(w[,i])),
                                    as.double(f@h), as.double(f@g),
                                    as.integer(f@L), as.double(Vout),
                                    PACKAGE="wavelets")[[7]])
                        }, w=DWj, v=DVj, f=filter, Vout=Vout)
          SVj <- sapply(1:n.series,
                        function(i,w,v,f,Vout){
                          return(.C("dwt_backward", as.double(w[,i]),
                                    as.double(v[,i]),
                                    as.integer(length(v[,i])),
                                    as.double(f@h), as.double(f@g),
                                    as.integer(f@L), as.double(Vout),
                                    PACKAGE="wavelets")[[7]])
                        }, w=SWj, v=SVj, f=filter, Vout=Vout)
        } else {
          Vout <- rep(0, length=dim(DWj)[1])
          DVj <- sapply(1:n.series,
                        function(i,w,v,f,k,Vout){
                          return(.C("modwt_backward", as.double(w[,i]),
                                    as.double(v[,i]), as.integer(k),
                                    as.integer(length(w[,i])),
                                    as.double(f@h), as.double(f@g),
                                    as.integer(f@L), as.double(Vout),
                                    PACKAGE="wavelets")[[8]])
                        }, w=DWj, v=DVj, f=filter, k=k, Vout=Vout)
          SVj <- sapply(1:n.series,
                        function(i,w,v,f,k,Vout){
                          return(.C("modwt_backward", as.double(w[,i]),
                                    as.double(v[,i]), as.integer(k),
                                    as.integer(length(v[,i])),
                                    as.double(f@h), as.double(f@g),
                                    as.integer(f@L), as.double(Vout),
                                    PACKAGE="wavelets")[[8]])
                        }, w=SWj, v=SVj, f=filter, k=k, Vout=Vout)
        }
      } else {
        if(method == "dwt"){
          DVj <- sapply(1:n.series,
                        function(i,w,v,f){
                          return(out <- dwt.backward(w[,i],v[,i],f))
                        }, w=DWj, v=DVj, f=filter)
          SVj <- sapply(1:n.series,
                        function(i,w,v,f){
                          return(out <- dwt.backward(w[,i],v[,i],f))
                        }, w=SWj, v=SVj, f=filter)
        } else {
          DVj <- sapply(1:n.series,
                        function(i,w,v,f,k){
                          return(out <- modwt.backward(w[,i],v[,i],f,k))
                        }, w=DWj, v=DVj, f=filter, k=k)
          SVj <- sapply(1:n.series,
                        function(i,w,v,f,k){
                          return(out <- modwt.backward(w[,i],v[,i],f,k))
                        }, w=SWj, v=SVj, f=filter, k=k)
        }
      }
      DWj <- matrix(rep(0, length=length(DVj)), ncol=n.series)
      SWj <- DWj
    }
    details[[j]] <- DVj
    smooths[[j]] <- SVj
  }
  
  # create mra object for output
  mra <- new("mra",
             D = details,
             S = smooths,
             filter = filter,
             level = as.integer(J),
             boundary = boundary,
             series = wt@series,
             class.X = wt@class.X,
             attr.X = wt@attr.X,
             method = method)

  return(mra)
}
