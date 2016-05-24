# modified from pbinomGC() in package tigerstats

binomGraphs  <- function(bound,region="below",size=100,prob=0.5,graph=TRUE,
                         xlab="x") {
  if (!is.numeric(bound)) stop("Specify one or two numerical boundaries")
  below <- grepl("^be[lf]",region,perl=TRUE)
  above <- grepl("^a[bf]",region,perl=TRUE)
  between <- grepl("^bet|^in",region,perl=TRUE)
  outside <- grepl("^out",region,perl=TRUE)
  if (length(bound)==1 & !(below | above)) stop("Specify region=\"below\" or
          region=\"above\"")
  if (length(bound)==2 & !(between | outside)) stop("Specify region=\"between\" or
          region=\"outside\"")
  if (length(bound)>2) stop("Specify one or two numerical boundaries")

  if (length(bound)==2 & bound[1]>bound[2])  bound <- rev(bound)

  sd <- sqrt(size*prob*(1-prob))
  
  if (below)  {
    area <- pbinom(bound,size=size,prob=prob)
    if (graph) {
    upper <- ceiling(max(qbinom(.9999,size=size,prob=prob),bound+0.1*sd))
    lower <- floor(min(qbinom(0.0001,size=size,prob=prob),bound-0.1*sd))
    nvals <- lower:upper
    Shading <- ifelse(nvals <= bound,"lightblue",NA)
    plot(nvals,dbinom(nvals,size=size,prob=prob),type="h",col=NA,axes=FALSE,
         xlab=xlab,ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
         main=paste("binom(",size,",",prob,") Distribution:\nShaded Area = ",round(area,3),sep=""))
    rect(nvals-0.5,rep(0,times=size+1),nvals+0.5,dbinom(nvals,size=size,prob=prob),
         col=Shading,border="black")
    axis(2)
    places <- c(lower,floor(bound),upper)
    axis(1,at=places,labels=c("",as.character(places[2]),""))
    }
  }

  if (above)  {
    area <- pbinom(bound,size=size,prob=prob,lower.tail=FALSE)
    if (graph) {
    upper <- ceiling(max(qbinom(.9999,size=size,prob=prob),bound+0.1*sd))
    lower <- floor(min(qbinom(0.0001,size=size,prob=prob),bound-0.1*sd))
    nvals <- lower:upper
    Shading <- ifelse(nvals > bound,"lightblue",NA)
    plot(nvals,dbinom(nvals,size=size,prob=prob),type="h",col=NA,axes=FALSE,
         xlab=xlab,ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
         main=paste("binom(",size,",",prob,") Distribution:\nShaded Area = ",round(area,3),sep=""))
    rect(nvals-0.5,rep(0,times=size+1),nvals+0.5,dbinom(nvals,size=size,prob=prob),
         col=Shading,border="black")
    axis(2)
    places <- c(lower,floor(bound)+1,upper)
    axis(1,at=places,labels=c("",as.character(places[2]),""))
    }
  }
  
  if (between)  {
    area <- pbinom(bound[2],size=size,prob=prob)-pbinom(bound[1]-1,size=size,prob=prob)
    if (graph) {
    upper <- ceiling(max(qbinom(.9999,size=size,prob=prob),bound+0.1*sd))
    lower <- floor(min(qbinom(0.0001,size=size,prob=prob),bound-0.1*sd))
    nvals <- lower:upper
    Shading <- ifelse((bound[1] <= nvals & nvals <= bound[2]),"lightblue",NA)
    plot(nvals,dbinom(nvals,size=size,prob=prob),type="h",col=NA,axes=FALSE,
         xlab=xlab,ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
         main=paste("binom(",size,",",prob,") Distribution:\nShaded Area = ",round(area,3),sep=""))
    rect(nvals-0.5,rep(0,times=size+1),nvals+0.5,dbinom(nvals,size=size,prob=prob),
         col=Shading,border="black")
    axis(2)
    places <- c(lower,floor(bound[1]),floor(bound[2]),upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  if (outside)  {
    area <- pbinom(bound[2],size=size,prob=prob,lower.tail=FALSE)+pbinom(bound[1]-1,size=size,prob=prob)
    if (graph) {
    upper <- ceiling(max(qbinom(.9999,size=size,prob=prob),bound+0.1*sd))
    lower <- floor(min(qbinom(0.0001,size=size,prob=prob),bound-0.1*sd))
    nvals <- lower:upper
    Shading <- ifelse(bound[1] <= nvals & nvals <= bound[2],NA,"lightblue")
    plot(nvals,dbinom(nvals,size=size,prob=prob),type="h",col=NA,axes=FALSE,
         xlab=xlab,ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
         main=paste("binom(",size,",",prob,") Distribution:\nShaded Area = ",round(area,3),sep=""))
    rect(nvals-0.5,rep(0,times=size+1),nvals+0.5,dbinom(nvals,size=size,prob=prob),
         col=Shading,border="black")
    axis(2)
    places <- c(lower,floor(bound[1])-1,floor(bound[2])+1,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }

}#end of binomGraphs