
#Function to plot the peeling/pasting trajectory and also identify other options
#for boxes that result from removing dimensions from more complicated boxes.

#(currently - 2008-10-22) has three options:  Plotting just the trajectory, 
#plotting "dominating" points that extent the tradeoff frontier, and plotting 
#contour lines of constant dimensionality.  

`trajplot` <-	function(infolist, coverage=TRUE, margtrajs, 
								colvect=c("red", "blue", "purple", "brown", "forest green"),
								trajplot_xlim, trajplot_ylim){

#infolist is object like that output from traj.info  - multiple lists, each of length trajectory  
#coverage indicates whether to plot coverage or support on x-axis
#margtrajs is ...
#colvect is the sequence of colors to use for points - they change when the
#dimensionality of the box changes - but the mapping is not one to one


#Turn off the warnings that will crop up from doing linear predictions with 
#only two points
currwarn <- getOption("warn")
options(warn=-1)

ltraj <- length(infolist$dimlist)

tdims <- pcvect <- vector(length=length(infolist$dimlist))

lcolvect <- length(colvect)

#tdims: how many dimensions are restricted for that box:
#pcvect: which color to use, as a function of restricted dimensions
for (i in 1:length(infolist$dimlist)){
  tdims[i] <- sum(infolist$dimlist[[i]]$either) #total number of dimensions rstr
  pcvect[i] <- colvect[tdims[i]%%lcolvect+1]    #color associated with that num
}

options("device")$device(width=7)

#Basic plots for either coverage or support:
if(coverage==TRUE){
  xax <- infolist$marcoverage
  plot(xax,infolist$y.mean, 
			xlab="Coverage", ylab="Density", 
			xlim=trajplot_xlim, ylim=trajplot_ylim,
			main="Peeling trajectory",
			col=pcvect, pch=19)
} else if(coverage==FALSE){
  xax <- infolist$mass
  plot(xax,infolist$y.mean, xlab="Support",ylab="Density", main="Peeling trajectory",col=pcvect)
} else (stop("coverage argument must be TRUE or FALSE"))

#This section should plot dimension contours
 bringToTop(-1)


#These statistics only guaranteed for coverage-oriented stats, 
#so we limit the option:
if(coverage){ 

	contours <- contourmkr(margtrajs) #contourmkr
	points(x=contours[[1]][1,1], y=contours[[1]][1,2], pch=15, cex=1, col="black")

  wht2plot <- readline(cat("Would you like to plot dimension contours, new dominating points,","\n",
    "or just continue on and pick boxes to inspect?","\n",
    "Enter 'dims','dom', or 'n')","\n"))
  
  if(wht2plot=="dims"){
  
    allpts <- matrix(ncol=4, nrow=0)
  
    for (i in 1:length(contours)){
  
      colind <- (i%%lcolvect)*(i!=lcolvect)+lcolvect*(i==lcolvect)
      
      points(contours[[i]],type="b",col=colvect[colind],pch=3, cex=.5)
      
      allpts <- rbind(allpts,cbind(contours[[i]],i-1,c(1:ltraj)))
      #we're cbinding the stats, plus the dimensionality + original box number
   
    }
    
  }
  
  
  ###RIGHT HERE SNIPPED OUT MORE ELEGANT CODE THAT DIDN'T QUITE WORK
  ###PASTED IN (COMMENTED OUT) AT END OF FUNCTION FOR FUTURE REFERENCE
  
  if(wht2plot=="dom"){  
  #plot only those points that extend the tradeoff frontier in cov-dens space
    
    if(length(contours)==1){
      cat("\n","Sorry, only one dimension, can't do dominated points.","\n")
      flush.console()
    } else{
    
    allpts <- matrix(ncol=4, nrow=0)
  
    for (i in 2:length(contours)){
  
      allpts <- rbind(allpts,cbind(contours[[i]],i-1,c(1:ltraj)))
      #we're cbinding the stats, plus the dimensionality + original box number
   
    }
  
  
    
    masbetters <- vector(length=nrow(allpts))  #oh... should end up being all true.
    
    for (i in 1:(length(xax)-1)){
      
      t1 <- (allpts[,1] > xax[i+1])
      t2 <- (allpts[,1] <=  xax[i])
      
      ininterval <- (t1 & t2) #The points that are in the interval  
      
      tdata <- data.frame(xax=xax[i:(i+1)],y=infolist$y.mean[i:(i+1)])
      tlm   <- lm(y~xax,tdata)
      
      preds <- predict.lm(tlm, data.frame(xax=allpts[,1]))
      betters <- (preds <= allpts[,2])
      betterses <- betters & ininterval 
      masbetters <- (masbetters | betterses)
      
    }
    
    points(allpts[masbetters,1:2],col=colvect[allpts[masbetters,3]%%5+1],pch=3,cex=.5)
  
    }
  
  }
}

cat("Now please select candidate boxes for inspection by clicking near them.","\n")
cat("After you have selected as many as you like, right click and select STOP.","\n","\n")


flush.console()


#Have to change what's passed to the identify() command dependent on what we 
#need to be able to identify

if(wht2plot=="dom"){

  labely <- c(c(1:ltraj),allpts[masbetters,4])
  idmat <- rbind(cbind(xax, infolist$y.mean),allpts[masbetters,1:2])
  idpts <- identify(idmat,labels=labely)
  idpts <- unique(labely[idpts])

} else if(wht2plot=="dims"){
  
  idmat <- rbind(cbind(xax, infolist$y.mean),allpts[,1:2])
  labely <- rep(c(1:ltraj),length(idmat))
  idpts <- identify(idmat,labels=labely)
  idpts <- unique(labely[idpts])
} else{
  idpts <- identify(xax,infolist$y.mean)
}

bringToTop(-1)
#To id points from removed dims: [NEED TO REMOVE OVERLAP]
#idpts <- identify(allpts[,1:2])
#
#for (i in idpts){
#
#  cind <- allpts[i,3]+1
#  
#  print(which(contours[[cind]]==allpts[i,1:2]))
#
#}

options(warn=currwarn) #return the warning setting to whatever it was.

return(idpts)



}




  #This section converts things into a [Maxd * nboxes] by 3 matrix with coordinates 
  #and dimension size in the third col
  
  #points(allpts[,1:2],col=colvect[allpts[,3]%%5+1],pch=3,cex=.5)
  
  #THIS SECTION IS DESIGNED TO FILTER AND FIND ONLY THOSE POINTS that
  #'DOMINATE' in the sense of extending the frontier, but with equal or 
  #lower dimensionality
  
  #THIS SECTION IMMEDIATELY BELOW ALMOST WORKS BUT NOT QUITE - SACRIFICING ELEGANCE FOR EVALUATION TIME, SKIP DOWN
  
  #masinterval <- vector(length=nrow(allpts))  #oh... should end up being all true.
  #
  #for (i in 1:length(xax)){
  #  
  #  t1 <- (allpts[,1] >= xax[i+1])
  #  t2 <- (allpts[,1] <  xax[i])
  #  
  #  ininterval <- (t1 & t2) #The points that are in the interval  
  #  
  #  masinterval <- (ininterval | masinterval)
  #  
  #  tdata <- data.frame(xax=xax[i:(i+1)],y=infolist$y.mean[i:(i+1)])
  #  tlm   <- lm(y~xax,tdata)
  #  
  #  if(any(ininterval)){
  #    
  #    preds <- predict.lm(tlm, data.frame(xax=allpts[ininterval,1]))
  #    betters <- (preds >= allpts[ininterval,2])
  #  
  ###BAD    masbetters <- (masbetters | betters)
  #  
  #  }
  #
  #}
  #
  ###BADpoints(allpts[masininterval,1:2][betters,1:2],col=colvect[allpts[,3]%%5+1],pch=3,cex=.5)
  