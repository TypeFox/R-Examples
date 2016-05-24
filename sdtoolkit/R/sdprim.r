
#This version corrects the argument passing problems used by in which 
#reproducibility statistics arguments were not making it into prim.traj

sdprim <-
function(x, y=NULL, thresh=NULL, peel.alpha=0.1, paste.alpha=0.05,
                 mass.min=0.001, pasting=TRUE, box.init=NULL,
                 coverage=TRUE, outfile="boxsum.txt", csvfile="primboxes.csv",
                 repro=TRUE,nbump=10,dfrac=.5, threshtype=">",
								 trajplot_xlim=c(0,1), trajplot_ylim=c(0,1),
								 peel_crit=1){
							
                 
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
#trajplot_[x/y]lim: each must be either a vector that can be passed as 
#              xlim/ylim arguments to plot, or NULL, which uses the 'old' version of 
#              plotting dim plots that doesn't necessarily include the full 
#  							boundaries of zero to 1 in each direction -- fills out the plot
#              	more but is slightly visually deceptive.
#peel_crit: peeling criteria to use -- 1 or 2.  1 is the default, 2 is the 
							  #alternate in eq 14.5 in Friedman and Fisher 1999, which 
								#a larger density increase if a larger box is removed
								#3 is minimizing mean of peeled box (eq 14.6).
				 
####GPL-3 notes:

cat("Welcome to the Scenario Discovery toolkit","\n","\n")

cat("Copyright (C) 2009  Evolving Logic","\n")
nicecat("This program comes with ABSOLUTELY NO WARRANTY; for details type 'sdwarranty()' when you are back at the R command prompt.  This is free software, and you are welcome to redistribute it under certain conditions; when you are not in this dialogue, type 'RShowDoc('COPYING')' for the complete GNU General Public License, which states these conditions.")

cat("\n")
  
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
		
	trajplot_xlimbad <- FALSE
	trajplot_ylimbad <- FALSE
# check format of trajplot_[x/y]lim and adjust if necessary:
	if(!is.null(trajplot_xlim)){
		if(is.numeric(trajplot_xlim)){
			if(length(trajplot_xlim)!=2){
				trajplot_xlimbad <- TRUE
			}
		}
	}
	# check format of trajplot_[x/y]lim and adjust if necessary:
	if(!is.null(trajplot_ylim)){
		if(is.numeric(trajplot_ylim)){
			if(length(trajplot_ylim)!=2){
				trajplot_ylimbad <- TRUE
			}
		}
	}
		
	if(trajplot_xlimbad){
		stop("trajplot_xlim has unacceptable argument: must be NULL or numeric vector
		     of length 2.")
	}
	
	if(trajplot_ylimbad){
		stop("trajplot_ylim has unacceptable argument: must be NULL or numeric vector
		     of length 2.")
	}
	
	xtemp <- x
	ytemp <- y
	d <- ncol(x)
	
	ranges	<- apply(xtemp,2,range)
	spans 	<- ranges[2,]-ranges[1,]
	
	## Onto actually doing stuff:
  keepgoing <- TRUE

  boxseq <- list()
  i = 0

  while (keepgoing){

    i <- i + 1	
  
#		#Drop dimensions with no variation:
	
#		ranges	<- apply(xtemp,2,range)
#		spans 	<- ranges[2,]-ranges[1,]
  
#		zero_spans <- which(spans==0)
#		if(length(zero_spans)>0){
#			cat("Found one or more columns with no variation, which will be dropped. They are: \n")
#			cat(colnames(xtemp)[zero_spans],"\n")
#			flush.console()
#			xtemp 			<- xtemp[,-zero_spans]
#			ranges  <- ranges[,-zero_spans]
#			spans   <- spans[-zero_spans]
#		}
		
		#Actually call peeling operations:
    tempbox <- prim.traj(xtemp, ytemp, box.init=box.init, peel.alpha=peel.alpha, 
      paste.alpha=paste.alpha, mass.min=mass.min, threshold=threshold, 
      pasting=pasting, verbose=verbose, threshold.type=threshold.type, 
      paste.all=paste.all,coverage=coverage,showbounds=showbounds,style=style,
      npts=npts,ninter=ninter, repro=repro, nbump=nbump, dfrac=dfrac,
			trajplot_xlim=trajplot_xlim, trajplot_ylim=trajplot_ylim,
			peel_crit=peel_crit)

    if(tempbox[1]!="done"){ #prim.traj will prompt to stop covering
    
      bxs <- tempbox$box
      rellow <- (bxs[1,]- ranges[1,])/spans
      relhi  <- (bxs[2,] - ranges[1,])/spans
    
      tempbox$relbox <- rbind(pmax(rellow,0.0),pmin(relhi,1.0))
      #The pmax's and mins are necessary because the other function
      #expands the initial box just to be safe...
    
      boxseq[[i]] <- tempbox

      #covering: brind identifies the points in the box, which are then removed
      #brind <- in.box(x=xtemp, box=tempbox$box, ncol(x), boolean=TRUE)
      brind <- in.box(x=xtemp, box=tempbox$box, ncol(xtemp), boolean=TRUE)
			xtemp <- xtemp[!brind,]
      ytemp <- ytemp[!brind]
      
      if(isTRUE(all.equal(ytemp,rep(0,length(ytemp))))){
     
        cat("No more nonzero points remaining, additional covering prohibited.","\n","\n")
        keepgoing <- FALSE
     
      } else {
  
        responsegood <- FALSE
    
        while(!responsegood){
    
        goon <- readline(cat("Continue covering? (\"y\" or \"n\")","\n"))
    
        if(goon!="y" & goon!="n"){
          cat("Please enter \"y\" or \"n\"","\n")
        } else {responsegood <- TRUE}
        
        }
        
        if(goon=="n"){keepgoing <-FALSE}

      }

    }
    
    else {keepgoing <- FALSE}

  }

  olm <- olaps(x,y,boxseq)

  estats <- uberstats(boxseq,npts,ninter,d)

  attr(boxseq,"estats") <- estats

  attr(boxseq,"olaps")  <- olm
  

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
  
  #autoprintout/writeout:
  
  seqinfo(boxseq,outfile)
  
  cat("\n")
  printout <- readline(cat("Would you like to print the boxes to csv format for export to CARs?","\n","(Enter \"y\" or \"n\")","\n"))
  
  
  if(printout=="y"){
  cat("The default filename to be printed is currently:",csvfile,"\n")
  filename <- readline(cat("If this is ok, enter 'y', otherwise enter the preferred filename (no quotes).","\n"))
  
  if(filename!="y"){
    csvfile <- filename
  }
  csvboxes(boxseq,outfile=csvfile)
  }
    
  cat("\n")
  return(boxseq)

}

