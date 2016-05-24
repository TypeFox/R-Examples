#' @title Chi Square Resampler (One at a Time)

#' @description An app to illustrate use of the chi-square statistic to test
#' for a relationship between two categorical variables.  The P-value is computed
#' by resampling, and the resamples are done one at a time.  A histogram of resampled chi-square statistics is displayed
#' after each resample, and summary information is output to the
#' console.
#' 
#' @rdname ChisqSimSlow
#' @usage ChisqSimSlow(form,data,effects=c("random","fixed"))
#' @param form a formula of the form ~x+y.  When using
#' fixed effects (see below for explanation), x should be the
#' variable that is considered the predictor variable.
#' @param data A data frame from which x and y are drawn.
#' @param effects When effects="fixed", the resampling is performed under
#' the conditon that the row sums in the resampled two-way table (with x
#' for rows) are the same as the row sums in the twoway table based on the
#' original data.  When effects="random", then both row and column sums
#' in the resampled table may vary: only the sum of the counts is
#' constant.  (Note:  in the resampling procedure for chisq.test
#' in the stats package of R, both row and column sums are
#' required to equal the corresponding sums for the original data.)
#' @return Graphical and numerical output
#' @export
#' @author Homer White \email{hwhite0@@georgetowncollege.edu}
#' @examples
#' \dontrun{
#' ChisqSimSlow(~weather+crowd.behavior,data=ledgejump,effects="fixed")
#' }
ChisqSimSlow <- function(form,data,effects=c("random","fixed"))  {
  
  #form is of the form ~var1+var2, both factors
  #best if var1 is explanatory, var2 is response
  #data is the data frame containing these variables
  #effects== "fnixed" performs resampling subject to constraint that row total
  #are fixed to observed row totals
  
  effects <- match.arg(effects)
  
  #Some utility functions
  rowsampler <- function(x,p)  {rmultinom(1,size=sum(x),prob=p)}  #used in slow sampling with fixed effects
  
  table.samp <- function(x)  {  #used in slow sampling with fixed effects
    nullprobs <- colSums(x)/sum(x)
    resamp <- t(apply(x,1,rowsampler,p=nullprobs))
    rownames(resamp)  <- rownames(x)
    colnames(resamp) <- colnames(x)
    as.table(resamp)
  }
  
  slowsamp <- function(x,effects)  {
    rnam <- rownames(x)
    cnam <- colnames(x)
    sofar <- matrix(0,nrow=nrow(x),ncol=ncol(x))
    rownames(sofar) <- rnam
    colnames(sofar) <- cnam
    callquick <- FALSE
    
    quickfin <- function(orig,unfin,effects)  {   #Quick finish for fixed effects
      
      if(effects=="fixed")  {
        first <- which(rowSums(unfin)<rowSums(orig))[1]
        rowleft <- rowSums(orig)[first]-rowSums(unfin)[first]
        for (j in 1:rowleft)  {
          item <- sampler(x)
          unfin[first,item] <- unfin[first,item]+1
        }
        
        if (first<nrow(x))  {
          for (i in (first+1):nrow(x))  {
            probs <- colSums(x)/sum(colSums(x))
            unfin[i,] <- rmultinom(1,size=sum(x[i,]),prob=probs)
          }
        }
        return(unfin)
      }
      
      if (effects=="random")  {
        numbleft <- sum(orig)-sum(unfin)
        return(unfin+rtabsamp(orig,numbleft))
      }
    }#end of quickfin
    
    
    sampler <- function(x)  {  #used in slow sampling when effects are fixed
      #x is a two-way table
      #sample at random from column-variable distribution given by table x, assume no relationship
      cdf.col <- cumsum(colSums(x))/sum(colSums(x))
      samp <- which(runif(1)<cdf.col)[1]
      samp
    }
    if (effects=="fixed")  {
      for (i in 1:nrow(x))  {
        if (callquick==TRUE)  break
        for (j in 1:rowSums(x)[i])  {
          w <- readline("Get an item, press Return.  To finish quickly, enter any character.")
          if (w=="")  {
            item <- sampler(x)
            cat("The sampled item is",rnam[i],cnam[item],"\n")
            sofar[i,item] <- sofar[i,item]+1
            cat("So far, the table looks like \n")
            print(sofar)
          } else {callquick <- TRUE;break}
        }
      }
      
      if (callquick==TRUE)  {
        sofar <- quickfin(orig=x,unfin=sofar,effects=effects)
      }
      
      return(sofar)
    }#end of effects==finxed
    
    if(effects=="random")  {
      for (i in 1:nrow(x))  {
        if (callquick==TRUE)  break
        for (j in 1:rowSums(x)[i])  {
          w <- readline("Get an item, press Return.  To finish quickly, enter any character.")
          if (w=="")  {
            item <- rtabsamp(x,1)  #retunrs table for sampling just one item
            rightrow <- which(rowSums(item)==1) #find row and column of th eitem
            rightcol <- which(colSums(item)==1)
            cat("The sampled item is",rnam[rightrow],cnam[rightcol],"\n")
            sofar <- sofar+item
            cat("So far, the table looks like \n")
            print(sofar)
          } else {callquick <- TRUE;break}
        }
      }
      
      if (callquick==TRUE)  {
        sofar <- quickfin(orig=x,unfin=sofar,effects=effects)
      }
      
      return(sofar)
    }#end of effects == random
  }#end of slowsamp
  
  exp.counts <- function(x)  (rowSums(x)%*%t(colSums(x)))/sum(x)
  
  #Resampler for "random effects":  row sums not fixed.
  rtabsamp <- function(x,n)  {
    #x is the table of observed counts
    #n is the number of items to resample
    #n=sum(x) amounts to one quick bootstrap resample
    
    expected <- exp.counts(x)
    probs <- expected/sum(x)
    resamp.tab <- rmultinom(1,size=n,prob=probs)
    resamp.tab <- matrix(resamp.tab,nrow=nrow(x))
    rownames(resamp.tab) <- rownames(x)
    colnames(resamp.tab) <- colnames(x)
    return(resamp.tab)
  }
  
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
  prsd <- with(data,ParseFormula(form))
  expname <- as.character(prsd$rhs)[2]
  respname <- as.character(prsd$rhs)[3]
  
  explanatory <- data[,expname]
  response <- data[,respname]
  
  x <- xtabs(~explanatory+response)
  
  N <- 50  #number of times player is liable to play
  overprob <- 0.025 #desired prob of going over xlimits and ylimits I will set
  q.needed <- 1+log(1-overprob)/N  #Using Poisson approx to binomial
  deg.freedom <- (nrow(x)-1)*(ncol(x)-1)
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
  expected <- exp.counts(x)
  rownames(expected)  <- rownames(x)
  colnames(expected) <- colnames(x)
  cat("If there is no relationship then the expected counts are:","\n")
  print(round(expected,1))
  stat <- sum((x-expected)^2/expected)
  cat("The observed value of the chisq statistic is:",round(stat,2),"\n")
  
  chisq.stats <- numeric()
  
  ymax.bar <- 1.3*max(x)
  barplot(t(x),beside=TRUE,
          legend.text=colnames(x),
          ylim=c(0,ymax.bar),xlab=expname,
          args.legend = list(x = "topleft",title=respname,cex=0.75),
          col=rainbow(ncol(x)),
          main=paste(expname,"and",respname))
  
  procedure <- readline("Enter s for slow sample, f for one fast sample, q to quit.")
  
  #Make sure they enter something:
  acceptable <- (procedure %in%c("s","f","q"))
  while(!acceptable)  {
    procedure <- readline("Hey!  I said:  Enter s for slow sample, f for one fast sample, q to quit.")
    acceptable <- (procedure %in%c("s","f","q"))
  }
  
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
    
    
    if (procedure=="s") {
      
      resamp <- slowsamp(x,effects=effects)
      expected <- exp.counts(resamp)
      newstat <- sum((resamp-expected)^2/expected)
      chisq.stats <- c(chisq.stats,newstat)
      
    }
    
    if (procedure=="f")  {
      if (effects=="fixed") resamp <- table.samp(x) else resamp <- rtabsamp(x,sum(x)) 
      expected <- exp.counts(resamp)
      newstat <- sum((resamp-expected)^2/expected)
      chisq.stats <- c(chisq.stats,newstat)
      
    }
    
    #output results to console:
    cat("\n")
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
    procedure <- readline("Enter s for slow sample, f for one fast sample, q to quit.")
    
    #Make sure they enter something:
    acceptable <- (procedure %in%c("s","f","q"))
    while(!acceptable)  {
      procedure <- readline("Hey!  I said:  Enter s for slow sample, f for one fast sample, q to quit.")
      acceptable <- (procedure %in%c("s","f","q"))
    }
  }
  
  
  
  return(cat("All done!"))
}#end of chisqSimSlow
