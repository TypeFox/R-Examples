plot.cumrates <-
function (x,ptrans,title,...)
{ # ptrans is a vector of character variables denoting the transitions selected, e.g. "NJ" 
   if (!inherits(x, "cumrates"))
        stop("'x' must be a 'cumrates' object")
  cumrates <- x
  if (missing(ptrans)) { print ("ptrans is missing"); return()}
  if (missing (title) ) title <- ""
  if (is.null(title)) title <= "" 	

  removed <- cumrates$D
  irate <- cumrates$irate
  cumh <- cumrates$predicted
  Lambda.oe <- cumrates$oeCum
  trans88 <- attr(removed,"param")$transitions$ODN
  plotrates <- which (trans88==ptrans)
# ========  Step 2: PLOT TRANSITION RATES ========
#           Plot cumulative hazard rates
#          Plot 2: cumulative hazard (as Putter)
if (length(ptrans)>0)   # no plot if ptrans = NULL
{ if (!exists("na")) 
	{ print ("No plot of Nelson-Aalen estimator since the computation of the estimator was skipped.")
	  return
	}
  colour <- c("red","darkgreen","blue","purple")
  ymax <- 0  
  nlegend <- length(ptrans)
  for (i in 1:attr(removed,"param")$ntrans)  ymax <- max(c(ymax,cumh[[i]]$na))
   namtransitions <- attr (removed,"param")$transitions$ODN
	if (irate %in% c(1,3))
  {  iz <- plotrates[1]  # first element of vector of transitions to be plotted
  	 x11 <- c(0,cumh[[iz]]$time)
   	 y11 <- c(0,cumh[[iz]]$na)
   	 y12 <- c(0,cumh[[iz]]$na-sqrt(cumh[[iz]]$var.aalen))
   	 y13 <- c(0,cumh[[iz]]$na+sqrt(cumh[[iz]]$var.aalen)) 
     plot(x11,y11,type="l",xlab="Age (years)",ylab="cumulative hazard",
        xlim=c(10,50), ylim=c(0,ymax),
        main=title,col=colour[1],axes=FALSE,lwd=2)
     lines (x11,y12,col=colour[1],lty=3,lwd=1)
     lines (x11,y13,col=colour[1],lty=3,lwd=1)
     axis (side=1,at=seq(10,50,by=5),labels=seq(10,50,by=5),cex.axis=0.8)
     axis (side=2,las=1,at=seq(0,ymax,by=0.5),
        labels=seq(0,ymax,by=0.5),cex.axis=0.8)
     # axis (side=2,las=1,at=seq(0,max(-log(sf1$surv)),by=0.5),labels=seq(0,max(-log(sf1$surv)),by=0.5),cex.axis=0.8)
     box()
     abline (h=seq(0,ymax+1,by=0.5),lty=2,col="lightgrey")
     abline (v=seq(10,50,by=5),lty=2,col="lightgrey")  # line at median age
     legend(10,ymax,namtransitions[plotrates],col=colour[1],
          lty=1,      # colour[1:nlegend],
          cex=0.9,bg="white",title="Ne-Aa estimator")
   
    if (length(ptrans)>1)  # more than one transition is plotted
     { 	for (ij in plotrates[2:length(plotrates)])
       { 
        x21 <- c(0,cumh[[ij]]$time)
   	    y21 <- c(0,cumh[[ij]]$na)
   	    y22 <- c(0,cumh[[ij]]$na-sqrt(cumh[[ij]]$var.aalen))
   	    y23<- c(0,cumh[[ij]]$na+sqrt(cumh[[ij]]$var.aalen))

        lines (x21,y21,col=colour[match(ij,plotrates)],lty=1,lwd=2)
        lines (x21,y22,col=colour[match(ij,plotrates)],lty=3)
        lines (x21,y23,col=colour[match(ij,plotrates)],lty=3)
        legend(10,ymax,namtransitions[plotrates],col=colour[1:nlegend],
          lty=1,      # colour[1:nlegend],
          cex=0.9,bg="white",title="Ne-Aa estimator")
       } } 
  }
    if (irate %in% c(2,3) & length(Lambda.oe)>1)
    { for (jj in 1:length(ptrans))
      { print (jj)
      	des <- as.numeric(as.character(attr(removed,"param")$transitions$DES[jj]))
      	or <-  as.numeric(as.character(attr(removed,"param")$transitions$OR[jj]))
      	lines (rownames(Lambda.oe),Lambda.oe[,des,or],col=colour[match(jj,plotrates)],lty=3,lwd=3)
      }
       legend (10,ymax/1.5,namtransitions[plotrates],col=colour[1:nlegend],
              lty=3,lwd=3,cex=0.9,bg="white",title="oe rates")
    }
 }

 return()
 }
