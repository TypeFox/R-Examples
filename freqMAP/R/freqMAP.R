`freqMAP` <-
function(dat,x,x.label,hw,cat.names=NULL,cat.short=NULL,num.samples=100000){

  if(is.null(cat.names)){
    #The list of unique cat names
    cat.names <- sort(unique(dat[,2]))
  }

  if(any(!(unique(dat[,2] %in% cat.names)))){
    stop("cat.names must contain at least all unique values of dat[,2]")
  }

  
  if(is.null(cat.short)){
    cat.short <- cat.names
  }else if(length(cat.names)!=length(cat.short)){
    stop("category shortforms must have same length as category names")
  }
  
  cat.ma <- cat.moving.average(dat=dat,x=x,hw=hw,cat.names=cat.names)
  names(cat.ma)[1] <- x.label

  #print(cat.ma)
  
  ll <- cat.freq.post(cat.ma=cat.ma,num.samples=num.samples)
  #print(ll)
  cat.ma <- ll$cat.ma
  post.samples <- ll$post.samples
  
  x <- list(cat.ma=cat.ma,post.samples=post.samples,
            cat.names=cat.names,cat.short=cat.short,
            hw=hw,x.label=x.label)
  class(x) <- c("freqMAP","list")
  x
}

