#' @title Graphical Calculator for Binomial Curve Probabilities

#' @description Shades desired areas under rectangles of probability histogram for binomial, 
#' returns numerical value of the area.
#' 
#' @rdname pbinomGC
#' @usage pbinomGC(bound,region="below",size=100,prob=0.5,graph=FALSE)
#' @param bound A numerical vector of length 1 or 2, range of shaded rectangles
#' @param region A character string.  Default is "below".  Possible values are "between" (when boundary consists of two numbers),
#' "below", "above", and "outside" (again when boundary consists of two numbers)
#' @param size Number of trials
#' @param prob Probability of success
#' @param graph produce graph?
#' @return Numerical value of probability.
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' #This gives P(X <= 6) for binom X with 10 trials, chance of success 0.70 on each trial:
#' pbinomGC(6,region="below",size=10,prob=0.70)
#' 
#' #This gives P(45 <= X <= 55), where X is binom with 100 trials,
#' #chance of success on each trial p = 0.50:
#' pbinomGC(c(45,55),region="between",size=100,prob=0.50)
#' 
#' #This gives P(X >= 7) = P(X > 6), for binom X with 10 trials,
#' #70% chance of success on each trial
#' pbinomGC(6,region="above",size=10,prob=0.7)
pbinomGC  <- function(bound,region="below",size=100,prob=0.5,graph=FALSE) {
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
         xlab="x",ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
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
         xlab="x",ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
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
         xlab="x",ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
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
         xlab="x",ylab="p(x)",xlim=c(lower-0.5,upper+0.5),
         main=paste("binom(",size,",",prob,") Distribution:\nShaded Area = ",round(area,3),sep=""))
    rect(nvals-0.5,rep(0,times=size+1),nvals+0.5,dbinom(nvals,size=size,prob=prob),
         col=Shading,border="black")
    axis(2)
    places <- c(lower,floor(bound[1])-1,floor(bound[2])+1,upper)
    axis(1,at=places,labels=c("",as.character(places[2:3]),""))
    }
  }
  
  return(area)

}#end of pbinomGC
