Sequences.ind.0 <-
function (d,namstates,absorb)
{  # d is object of class 'numeric' or 'Date'
	if (missing (absorb)) absorb <- NULL
	if (missing(namstates)) namstates <- NULL
	if (is.null(dim(d))) d <- matrix(d,nrow=1,ncol=length(d)) # if d is vector
    if (is.null(namstates)) {namstates=LETTERS[1:ncol(d)]
   	                        colnames(d) <- namstates}
    if(is.null(colnames(d))) colnames(d) <- namstates	

	if (!is.null(absorb)) {if(FALSE %in% (absorb %in% namstates)) absorb <- NULL
                           nn <- which (colnames(d) %in% absorb)
                           d[,nn] <-d[,nn] + sign(d[,nn]) * 0.1
                           for (i in 1:nsample)
                             {  for (jk in 1:length (nn))
                             	{  if (!is.na(d[i,nn[jk]])) d[i,] <- ifelse (d[i,]>d[i,nn[jk]],NA,d[i,])  }
                                d[,nn] <- trunc(d[,nn])
                             }
                           }
	first <- namstates[1]
	if (!is.character(first)) warning ("The first state is not a character.")
   if (is.vector(d)) d <- t(as.matrix(d)) # single observation (nsample=1)
   nsample <- nrow(d)
   d_s <- array(NA,dim=c(nrow(d),ncol(d)))
   sequence <- data.frame(d_s)
   path <- vector(mode="character",length=nsample)
   path <- rep(first,nsample)
   for (i in 1:nsample)  
     {  if ("FALSE"%in%is.na(d[i,]))  # vector has real values in addition to NA
     	 { zx <- sort (d[i,],decreasing=FALSE,na.last=TRUE) # needs to be data frame
     	   xx <- sum(!is.na(zx))     #  length(na.omit(zx))
    	   d_s[i,] <- as.numeric(zx)
    	   sequence[i,1:xx] <- names(zx)[1:xx]
    	   for (j in 1:xx) path[i] <- paste(path[i],names(zx)[j],sep="") } else
    	 { d_s[i,] <- as.numeric(d[i,]) }
     }
   #d_s <- data.frame(d_s)
   return (list(namstates=namstates,
                d=d_s,   # numeric
                path=path))
}
