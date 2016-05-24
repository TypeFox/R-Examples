`pavasmacof` <-
function(x, w = rep(1,length(x)), block=weighted.mean)
{

#---------------------- subroutines ---------------------------
  merge.block.up <- function(blocklist, blockvalues, x, w, i, block)
  {
    n <- length(blockvalues)
    nn <- 1:n
    ii <- which(i+1!=nn)
    blocklist[i,] <- c(blocklist[i,1],blocklist[i+1,2])
    indi <- blocklist[i,1]:blocklist[i+1,2]
    blockvalues[i] <- block(x[indi],w[indi])
    blocklist <- blocklist[ii,]
    if (length(ii) == 1) dim(blocklist)<-c(1,2)
    blockvalues <- blockvalues[ii]
    list(v = blockvalues, l = blocklist)
  }

  put.back <- function(n, blocklist, blockvalues)
  {
    x <- rep(0,n)
    nb <- length(blockvalues)
    for (i in 1:nb) x[blocklist[i,1]:blocklist[i,2]] <- blockvalues[i]
    return(x)
  }

  is.up.satisfied <- function(x,i) (i == length(x))||(x[i]<=x[i+1])
  is.down.satisfied <- function(x,i) (i == 1)||(x[i-1]<=x[i])

  weighted.median <- function(x,w=rep(1,length(x)))
  {
    ox <- order(x)
    x <- x[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0,w))
    up<-sum(w)-low
    df<-low-up
    repeat
    {
      if (df[k] < 0) k<-k+1
	else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
	       else return(x[k-1])
    }
  }

  weighted.pth.fractile <- function(x,w=rep(1,length(x)),a=1,b=1)
  {
    ox <- order(x)
    x <- x[ox]
    w <- w[ox]
    k <- 1
    low <- cumsum(c(0,w))
    up <- sum(w)-low
    df <- a*low-b*up
    repeat
    {
	if (df[k] < 0) k<-k+1
	  else if (df[k] == 0) return((w[k]*x[k]+w[k-1]*x[k-1])/(w[k]+w[k-1]))
		 else return(x[k-1])
    }
  }





#-------------------- end subroutines -------------------------

nblock <- n <-length(x) 
blocklist<-array(1:n,c(n,2)) 
blockvalues<-x
active<-1
repeat{
	if (!is.up.satisfied(blockvalues,active)) {
		blockmerge<-merge.block.up(blocklist,blockvalues,x,w,active,block)
		blockvalues<-blockmerge$v; blocklist<-blockmerge$l
		nblock<-nblock-1
		while (!is.down.satisfied(blockvalues,active)) {
			blockmerge<-merge.block.up(blocklist,blockvalues,x,w,active-1,block)
			blockvalues<-blockmerge$v; blocklist<-blockmerge$l; 
			nblock<-nblock-1; active<-active-1;
			}
		}
	else if (active == nblock) break() else active<-active+1
	}	
put.back(n,blocklist,blockvalues)
}

