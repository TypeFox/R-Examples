#' @title Chi Square Resampler (One at a Time) for Goodness-of-Fit

#' @description An app to illustrate use of the chi-square statistic to test
#' for goodness of fit.  The P-value is computed
#' by resampling, and the resamples are done one at a time.  A histogram of resampled chi-square statistics is displayed
#' after each resample, and summary information is output to the
#' console.
#' 
#' @rdname SlowGoodness
#' @usage SlowGoodness(x,p)
#' @param x a one-dimensional table, or a vector of observed counts
#' @param p vector of null probabilities
#' @return Graphical and numerical output
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' throws <- c(one=8,two=18,three=11,four=7,five=9,six=7)
#' SlowGoodness(throws,p=rep(1/6,6))
#' }
SlowGoodness <- function(x,p)  {
  
  
  #Some utility functions
   
  BinFinder <- function(val,breaks)  {
    #breaks a numeric vector of strictly increasing values
    if (val < breaks[1]) return(c(-Inf,breaks[1]))
    if (val >=breaks[length(breaks)]) return(c(breaks[length(breaks)],Inf))
    left <- max(which(breaks<=val))
    right <- min(which(val<breaks))
    return(c(breaks[left],breaks[right])) 
  }
  
  RightRow <- function(counts,val)  {
    n <- nrow(counts)
    breaks <- c(counts[,1],counts[n,2])
    rightrow <- which(counts[,1]==BinFinder(val,breaks)[1])
    return(rightrow)
  }
  
  UpdateCounts <- function(counts,new)  {
    #counts is a matrix with 3 columns:
    #LH sides of bins, RH sides of bins, count in bin
    #new is the new value
    rightrow <- RightRow(counts,new)
    counts[rightrow,3] <- counts[rightrow,3]+1
    return(counts)
  }
  
  UpdateHist <- function(counts,prev,new,oldcolor,newcolor)  {
    #make the bin of prev all oldcolor
    if(length(chisq.stats)>1) {
      rr.prev <- RightRow(counts,prev)
      rect(xleft=counts[rr.prev,1],
           ybottom=0,
           xright=counts[rr.prev,2],
           ytop=counts[rr.prev,3],
           col=oldcolor,border="black")
    }
    
    #add a distinct rectangle to top of bin for new value
    rr.new <- RightRow(counts,new)
    rect(xleft=counts[rr.new,1],
         ybottom=counts[rr.new,3],
         xright=counts[rr.new,2],
         ytop=counts[rr.new,3]+1,
         col=newcolor,border="black")
    
    #draw vertical line at stat just in case it got obscured:
    lines(x=c(stat,stat),y=c(0,ymax),lwd=2)
  }
  
  ResetHist <- function(hist,new,oldcolor)  {
    oldxmax <- max(hist[,2])
    newxmax <- floor(new)+1
    xmax <<- newxmax  #write upstairs just in case it's needed later
    newcountcol <- c(hist[,3],rep(0,newxmax-oldxmax))
    newleftcol <- 0:(newxmax-1)
    newrightcol <- 1:newxmax
    hist.info <- cbind(newleftcol,newrightcol,newcountcol)
    
    #Set up the plot:
    plot(0,0,col=rgb(0,0,1,0),xlim=c(0,newxmax),ylim=c(0,ymax),
         xlab="Resampled Chisq values",
         ylab="Count",
         main="Distribution of Resampled Chisq Values So Far")
    
    #add the rectangles, all oldcolor
    for(i in 1:nrow(hist.info))  {
      rect(xleft=hist.info[i,1],
           ybottom=0,
           xright=hist.info[i,2],
           ytop=hist.info[i,3],
           col=oldcolor,border="black")
    }
    
    return(hist.info)  
  }
  
  #___________________________
  #Ok, get started on processing the input
  
  x <- as.table(x)
  expected <- sum(x)*p
  
  N <- 50  #likely number of resamples before user tires
  overprob <- 0.025 #desired prob of going over xlimits and ylimits I will set
  q.needed <- 1+log(1-overprob)/N  #Using Poisson approx to binomial
  deg.freedom <- nrow(x)-1
  xmax <- ceiling(qchisq(q.needed,df=deg.freedom))  #set upper limit on horiz axis
  
  #Now go for upper limit on vertical axis (counts):
  left <- 0:(xmax-1)
  right <- 1:xmax
  leftprobs <- pchisq(left,df=deg.freedom)
  rightprobs <- pchisq(right,df=deg.freedom)
  binprobs <- rightprobs-leftprobs
  biggie <- max(binprobs)
  ymax <- ceiling(qbinom(1-overprob,size=N,prob=biggie))
  
  #Set up initial info for histogram:
  breaks <- 0:xmax
  m <- length(breaks)
  hist.info <- cbind(breaks[-m],breaks[-1],rep(0,(m-1)))
  prev <- 0
  new <- 0
  oldcolor <- "lightgreen"
  newcolor <- "red"
  
  #lay it out for the user:
  cat("The observed cell counts are:","\n")
  print(x)
  cat("\n")
  #now expected cell counts:
  expected <- sum(x)*p
  names(expected)  <- names(x)
  cat("If the Null is right then the expected counts are:","\n")
  print(round(expected,1))
  stat <- sum((x-expected)^2/expected)
  cat("The observed value of the chisq statistic is:",round(stat,2),"\n")
  
  chisq.stats <- numeric()
  
  ymax.bar <- 1.3*max(x)
  barplot(t(x),beside=TRUE,
          ylim=c(0,ymax.bar),xlab="", ylab="Frequency",
          main="Tally of Variable")
  
  procedure <- readline("Press Enter for one fast sample, enter q to quit.")
  
  #Make sure they enter something:
  
  if (procedure=="q")  return(cat("All Done!"))
  else {
    plot(0,0,col=rgb(0,0,1,0),xlim=c(0,xmax),ylim=c(0,ymax),
         xlab="Resampled Chisq values",
         ylab="Count",
         main="Distribution of Resampled Chisq Values So Far")
    
  }
  
  while(procedure != "q")  {
    
    max.height <- max(hist.info[,3])
    if(max.height >=ymax)  {#This person is persistant.  Give her more head room.
      cat("My, my, you ARE the persistent type!\n")
      cat("I'll raise the vertical axis so you can keep going.\n")
      ymax <- floor(1.5*ymax)      
      #Trick here:  reuse the ResetHist function
      ResetHist(hist.info,xmax-1,oldcolor)
    }
    
      resamp <-  t(rmultinom(1, size = sum(x), prob = p))
      colnames(resamp) <- names(x)
      newstat <- sum((resamp-expected)^2/expected)
      chisq.stats <- c(chisq.stats,newstat)
    
    #output results to console:
    cat("\n")
    resamp <- as.table(resamp)
    rownames(resamp) <- ""
    print(resamp)
    cat("For this resample the chisq statistic is",round(newstat,2),".")
    cat("Total resamples so far is",length(chisq.stats),"\n")
    perc <- 100*length(chisq.stats[chisq.stats >= stat])/length(chisq.stats)
    cat("\n")
    cat("Percentage of times the resampled statistics have exceeded observed chisq statistic",
        round(stat,2),"is",round(perc,1),"%")
    
    #now update the histogram
    #first deal with possibility that newstat exceeds old bounds:
    if(newstat>xmax)  hist.info <- ResetHist(hist.info,newstat,oldcolor)
    #Then update:
    new <- newstat 
    UpdateHist(hist.info,prev,new,oldcolor,newcolor)
    hist.info <- UpdateCounts(hist.info,new)
    prev <- new
    
    #Query for another resample:
    procedure <- readline("Press Enter for one fast sample, enter q to quit.")
    
  } # end of while procedure != q
  
  
  return(cat("All done!"))
  
}#end of SlowGoodness
