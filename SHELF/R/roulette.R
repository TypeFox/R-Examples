#' Elicit one set of probabilities using the roulette method.
#' 
#' Produces a graphics window with the roulette grid. The user clicks in the
#' window to allocate 'chips' to 'bins'. The elicited probability inside each
#' bin is the proportion of chips in each bin.
#' 
#' 
#' @param lower The lower limit on the x-axis of the roulette grid.
#' @param upper The upper limit on the x-axis of the roulette grid.
#' @param gridheight The maximum number of chips that can be allocated to a
#' single bin.
#' @param nbins The number of equally sized bins drawn between \code{Lo} and
#' \code{Up}.
#' @param round.end If set to \code{T}, empty bins and the uppermost non-empty
#' bin will be ignored. For example, with 20 chips in total, if the uppermost
#' non-empty bin is [70,80] and contains 1 chip, setting \code{round.end = F}
#' will result in an elicited probability P(X>80)=0, but setting
#' \code{round.end = T} will remove this judgement, instead only having
#' P(X>70)=0.05.
#' @return A list, with outputs 
#' \item{v }{ upper limits of
#' each bin.}
#' \item{p }{ cumulative probabilities for each
#' upper bin limit.}
#' @author Jeremy Oakley <j.oakley@@sheffield.ac.uk>
#' @examples
#' 
#' \dontrun{
#' x <- roulette()
#' # Then allocate chips to bins
#' 
#' # To fit distributions and see the results
#' myfit <- fitdist(vals = x$v, probs = x$p)
#' plotfit(myfit)
#' }
#' @export
roulette<-function(lower=0, upper=100, gridheight=10, nbins=10, round.end = F){
  
  chips<-rep(0,nbins)
  bin.width<-(upper-lower)/nbins
  bin.left<-seq(from=lower,to=upper-bin.width,length=nbins)
  bin.right<-seq(from=lower+bin.width,to=upper,length=nbins)
  xy<-list(x=lower,y=0)

  while(xy$x>= lower & xy$x<= upper){

      par(mfrow=c(1,1))
      plot(c(lower,upper),c(0,0),xlim=c(lower,upper),ylim=c(-1,max(gridheight,max(chips)+1)),type="l",ylab="",xaxp=c(lower,upper,nbins), main = paste("Total chips:", sum(chips)), xlab="Click outside the x-axis range to finish")
      for(i in 1:nbins){
        lines(c(bin.left[i],bin.left[i]),c(0,max(gridheight,max(chips)+1)),lty=3,col=8)
      }
      lines(c(bin.right[nbins],bin.right[nbins]),c(0,max(gridheight,max(chips)+1)),lty=3,col=8)

    for(i in 1:gridheight){
      lines(c(lower,upper),c(i,i), lty=3,col=8)
    }
      
    for(i in 1:nbins){
        if(chips[i]>0){
          rect(rep(bin.left[i],chips[i]),c(0:(chips[i]-1)),rep(bin.right[i],chips[i]),c(1:chips[i]),col=2)
        }
      }

    xy<-locator(n=1)
    if(length(xy)==0){xy=list(x=-1,y=-1)}
  
    if(xy$x > lower & xy$x <upper & xy$y < gridheight){
      index<-ceiling(xy$x/upper*nbins)
      chips[index]<-ceiling(max(xy$y,0))}
    }

  outputs<-list(v=bin.right, p=cumsum(chips)/sum(chips))
  if(round.end == T){
    index <- outputs$p>0 & outputs$p<1
    outputs <- list(outputs$v[index], outputs$p[index])
  }
  
  # Redraw final plot without x-axis label
  
  plot(c(lower,upper),c(0,0),xlim=c(lower,upper),ylim=c(-1,max(gridheight,max(chips)+1)),type="l",ylab="",xaxp=c(lower,upper,nbins), main = paste("Total chips:", sum(chips)), xlab="")
  for(i in 1:nbins){
    lines(c(bin.left[i],bin.left[i]),c(0,max(gridheight,max(chips)+1)),lty=3,col=8)
  }
  lines(c(bin.right[nbins],bin.right[nbins]),c(0,max(gridheight,max(chips)+1)),lty=3,col=8)
  
  for(i in 1:gridheight){
    lines(c(lower,upper),c(i,i), lty=3,col=8)
  }
  
  for(i in 1:nbins){
    if(chips[i]>0){
      rect(rep(bin.left[i],chips[i]),c(0:(chips[i]-1)),rep(bin.right[i],chips[i]),c(1:chips[i]),col=2)
    }
  }
  
  outputs
}    
