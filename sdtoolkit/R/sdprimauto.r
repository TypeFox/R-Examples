#A version of sdprim that will run on it's own, automatically selecting
#highest density 100% coverage box and highest coverage 100% density box
#At present, it will just stop after doing one box


sdprimauto <-
function(x, y=NULL, thresh=NULL, peel.alpha=0.1, paste.alpha=0.05,
                 mass.min=0.001, pasting=TRUE, box.init=NULL,
                 coverage=TRUE, outfile="boxsum.txt", csvfile="primboxes.csv",
                 repro=FALSE, nbump=10,dfrac=.5, threshtype=">", inc_dset=FALSE){
                 
#Arguments:
#x: input data
#y: output data (thresholded)
#box.init: IGNORE
#peel.alpha: peeling patience parameter
#paste.alpha: pasting parameter
#mass.min:  Minimum support to allow continued peeling of a box
#threshold: don't mess with this for now
#pasting: is currently automatic, regardless of what you put here (I think)
#verbose: old implementation prints out when done peeling and pasting
#threshold.type: don't mess with this for now
#coverage: logical, show plots and stats in terms of  coverage, rather than support
#outfile: file to print a copy of the console box summary information to - not a csv file, more easily readable
#csvfile: file to print csv box information for reading into CARs
#repro: logical, whether or not to produce reproducibility statistics (which are achieved by rerunning PRIM on resamplings of the dataset
#nbump:  given that repro==TRUE, how many resamplings should be used?
#dfrac:  given that repro==TRUE, what fraction of the data should be resampled each time?
         
####GPL-3 notes:


#
#cat("Welcome to the Scenario Discovery toolkit","\n","\n")
#
#cat("Copyright (C) 2009  Evolving Logic","\n")
#nicecat("This program comes with ABSOLUTELY NO WARRANTY; for details type 'sdwarranty()' when you are back at the R command prompt.  This is free software, and you are welcome to redistribute it under certain conditions; when you are not in this dialogue, type 'RShowDoc('COPYING')' for the complete GNU General Public License, which states these conditions.")
#
#cat("\n")
#  
#arguments which get passed on, but which should really just be removed, at least for now:
showbounds <- TRUE
style <- "ineq"
paste.all <- TRUE
threshold.type <- 1
verbose <- FALSE
threshold <- 0
  
##VARIOUS CHECKS:


  sdprimcall <- match.call()
  if(is.null(sdprimcall$y)){ #then y not specified, and remove the last column of x
    y <- x[,ncol(x)]
    x <- x[,1:(ncol(x)-1)]

  }

  #Threshold the data if not thresholded
  if(!is.null(thresh)){
    yreal <- y  #store the real valued y, just in case
    relop <- get(threshtype)
    y <- 1*(relop(y,thresh))
  }
  
  
  if(is.null(dim(x))){stop("Sorry, can't handle single column x matrix at this point.  But you shouldn't really need PRIM anyway.")} 

  npts <- nrow(x)
  ninter <- sum(y)

#  mass.min=ceiling(1/peel.alpha)

  #check there are nonzero y's
  if(isTRUE(all.equal(y,rep(0,npts)))){stop("Need at least some nonzero values in output vector y.")}
    
  keepgoing <- TRUE

  xtemp <- x
  ytemp <- y
  d <- ncol(x)
  
  ranges <- apply(x,2,range)
  spans <- ranges[2,]-ranges[1,]
  
  boxseq <- list()
  i = 0

  while (keepgoing){

    i <- i + 1

    tempbox <- prim.trajauto(xtemp, ytemp, box.init=box.init, peel.alpha=peel.alpha, paste.alpha=paste.alpha, mass.min=mass.min, threshold=threshold, pasting=pasting, verbose=verbose, threshold.type=threshold.type, paste.all=paste.all,coverage=coverage,showbounds=showbounds,style=style,npts=npts,ninter=ninter, repro=repro)


    if(tempbox[1]!="done"){ #prim.traj will prompt to stop covering
    
      bxs <- tempbox$box
      rellow <- (bxs[1,]- ranges[1,])/spans
      relhi  <- (bxs[2,] - ranges[1,])/spans
    
      tempbox$relbox <- rbind(pmax(rellow,0.0),pmin(relhi,1.0))
      #The pmax's and mins are necessary because the other function
      #expands the initial box just to be safe...
    
      boxseq[[i]] <- tempbox

      keepgoing <- FALSE #(autoaddition)

#      #cover
#  
#      brind <- in.box(x=xtemp, box=tempbox$box,ncol(x),boolean=TRUE)
#  
#      xtemp <- xtemp[!brind,]
#      ytemp <- ytemp[!brind]
#  
#      
#      if(isTRUE(all.equal(ytemp,rep(0,length(ytemp))))){
#     
#        cat("No more nonzero points remaining, additional covering prohibited.","\n","\n")
#        keepgoing <- FALSE
#     
      } #else {
#  
#        responsegood <- FALSE
#    
#        while(!responsegood){
#    
#        goon <- readline(cat("Continue covering? (\"y\" or \"n\")","\n"))
#    
#        if(goon!="y" & goon!="n"){
#          cat("Please enter \"y\" or \"n\"","\n")
#        } else {responsegood <- TRUE}
#        
#        }
#        
#        if(goon=="n"){keepgoing <-FALSE}
#
#      }
#
#    }
#    
#    else {keepgoing <- FALSE}

  }



  olm <- olaps(x,y,boxseq)

  estats <- uberstats(boxseq,npts,ninter,d)

  attr(boxseq,"estats") <- estats

  attr(boxseq,"olaps")  <- olm
  
  if(inc_dset==TRUE){  

    indmat <- matrix(nrow=nrow(x),ncol=length(boxseq))
    colnames(indmat) <- c(1:ncol(indmat))
  
    for (i in 1:length(boxseq)){
  
      indmat[,i] <- in.box(x,boxseq[[i]]$box,ncol(x),boolean=TRUE) #OR WHATEVER
      colnames(indmat)[i] <- paste("in.B",i,sep="")
  
    }
  
    inseq <- apply(indmat,1,any)
  
  #convert to 0-1 for dataset output
    indmat <- 1*indmat 
    inseq  <- 1*inseq
  
  
  
    augmented <- cbind(x,y,inseq,indmat)
    colnames(augmented) <- c(colnames(x),"output","in.seq",colnames(indmat))
  
    attr(boxseq,"data") <- augmented
  }
  
  #autoprintout/writeout:
  
#  seqinfo(boxseq,outfile)
  
#  cat("\n")
#  printout <- readline(cat("Would you like to print the boxes to csv format for# export to CARs?","\n","(Enter \"y\" or \"n\")","\n"))
  
  
#  if(printout=="y"){
#  cat("The default filename to be printed is currently:",csvfile,"\n")
#  filename <- readline(cat("If this is ok, enter 'y', otherwise enter the preferred filename (no quotes).","\n"))
  
#  if(filename!="y"){
#    csvfile <- filename
#  }
#  csvboxes(boxseq,outfile=csvfile)
#  }
    
#  cat("\n")
  return(boxseq)

}

