################Ordered FDR procedure
fdrOrd <- function(p,q=.01,ord=NULL,GD=FALSE){ #tune=1,,
  m <- length(p)
  alphaprime=q
  
  # sort by ord
  if(!is.null(ord))    o <- order(ord,decreasing=T)
  else { o <- 1:length(p) ; ord= o[length(p):1]}
  
  if(GD) alphaprime=alphaprime/(sum(1/((1+(1:length(p)))/(2-alphaprime))))
  else alphaprime=q
  
  u <- cumsum(p[o]>alphaprime)
  h <- rep(0,m)
  names(h) <- names(p)
  w=min(c(which((u*alphaprime/(1-alphaprime))/((1:m)-u)>q)-1,m))
  #w=min(c(which((u*tune*alphaprime/(1-alphaprime))/((1:m)-u)>q)-1,m))
  h[o[1:w]] <- 1
  h[p>alphaprime] <- 0

  out <- new("someMTP.object")  
  out @rej = h == 1
  out @p = p
  out @ord = ord
  out @idOrd= o
  out @MTP = "fdrOrd"
  out @GD = GD
  out @q = q
  out @k = NULL
  out @alpha = NULL
  out @alphaprime = alphaprime  
  
  return(out)
}

###############FDR by Benjamini Hochberg
# fdrBH <- function(p,q=.01,disp=TRUE) {
  # m <- length(p)
  # ps <- sort(p)
  # h <- rep(0,m)
  # ti <- 0
  # u <- ps<(q*(1:m)/m)
  # if(any(u))     ti <- ps[max(which(u))]
  # h[which(p<=ti)] <- 1
  # if(disp) {
    # cat(paste("Benjamini Hochberg FDR procedure\n ",length(p)," tests, q=",q,"\n threshold=",ti,"\n ",sum(h)," rejections\n\n",sep=""))
  # }
  # return(h==1)
# }