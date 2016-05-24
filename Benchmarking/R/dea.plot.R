# $Id: dea.plot.R 140 2015-05-15 21:48:02Z B002961 $
"dea.plot" <-
function(x, y, RTS="vrs", ORIENTATION="in-out", txt=NULL, add=FALSE, 
            wx=NULL, wy=NULL, TRANSPOSE = FALSE, fex=1, GRID=FALSE,
            RANGE=FALSE, param=NULL, ..., xlim, ylim, xlab, ylab)
# x er input 1 og y er iput 2 eller output.
# Hvis der flere varer i de to input/output bliver de lagt sammen som
# vaegtet sum med vaegte wx og wy; default vaegte som vaere simpel addition.
#
{
   rts <- c("fdh","vrs","drs","crs","irs","irs2","add","fdh+")
   if ( is.numeric(RTS) )  {
      cat("Number '",RTS,sep="")
      RTStemp <- rts[1+RTS] # the first fdh is number 0
      RTS <- RTStemp
      cat("' is '",RTS,"'\n",sep="")
   }
   RTS <- tolower(RTS)
   if ( !(RTS %in% rts) ) stop(paste("Unknown value for RTS:",RTS),quote=F)

   orientation <- c("in-out","in","out","graph")
   if ( is.numeric(ORIENTATION) )  {
      ORIENTATION_ <- orientation[ORIENTATION+1]  # "in-out" er nr. 0
      ORIENTATION <- ORIENTATION_
   }
   ORIENTATION <- tolower(ORIENTATION)
   if ( !(ORIENTATION %in% orientation) ) {
      stop(paste("Unknown value for ORIENTATION:",ORIENTATION))
   }

   if ( RTS=="fdh+" && ORIENTATION!="in-out" )
      stop("RTS=\"fdh+\" only works for ORIENTATION=\"in-out\"")

   if ( RTS=="add" && ORIENTATION!="in-out" )
      stop("RTS=\"add\" only works for ORIENTATION=\"in-out\"")

   if (TRANSPOSE) {
      x <- t(x)
      y <- t(y)
      if ( !is.null(wx) )  {
          if ( is.matrix(wx) )  { wx <- t(wx) }
      }
      if ( !is.null(wy) )  {
         if ( is.matrix(wy) )  { wy <- t(wy) }
      }
   }

   if ( is.matrix(x) && dim(x)[2] > 1 )  {  
         # x skal aggregeres
         if ( is.null(wx) ) { wx <- matrix(1, nrow=dim(x)[2] ,ncol=1) }
         x <- x %*% wx
   }
   if ( is.matrix(y) && dim(y)[2] > 1 )  {  
         # y skal aggregeres
         if ( is.null(wy) ) { wy <- matrix(1, nrow=dim(y)[2] ,ncol=1) }
         y <- y %*% wy
   }

   if ( add == FALSE ) {
      dots = list(...)
	   if (RANGE)  {
		   xlim <- 1.2*range(c(0,x)) +c(0,.01)
		   ylim <- 1.2*range(c(0,y)) +c(0,.01)
	   }
      if ( missing(xlim) ) xlim=c(0,1.2*max(x+.001))
      if ( missing(ylim) ) ylim=c(0,1.2*max(y+.0011))
      if ( missing(xlab) ) {xlab <- switch(ORIENTATION, 
                        "in"="x1", "out"="y1", "in-out"="X")}
      if (missing(ylab) ) {ylab <- switch(ORIENTATION, 
                        "in"="x2", "out"="y2", "in-out"="Y")}

      if ( RTS=="fdh+" ) {
         if ( is.null(param) )  {
            delta <- .15
            low <- 1-delta
            high <- 1+delta
         } else {
            if ( length(param) == 1 )  {
               low <- 1-param
               high <- 1+param
            } else {
               low <- param[1]
               high <- param[2]
            }
         }
         xlim <- c(low,high)*xlim
         ylim <- c(low,high)*ylim
      }

      # plot points with axes
      plot(x,y,xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",xlab=xlab,ylab=ylab,
                frame=FALSE,...)
      if ( GRID )  {
         grid(col="darkgray")
         box(col="grey")
      }
      if ( class(txt)=="logical" && txt )  {
         if ( class(x)=="matrix" )  {
            if ( !is.null(rownames(x)) )  {
               txt <- rownames(x)
            } else if ( !is.null(rownames(y)) )  {
               txt <- rownames(y)
            } else {
               txt <- 1:dim(x)[1]
            }
         } else {
            txt <- 1:length(x)
         }
      }
      if ( class(txt)!="logical" && length(txt) > 0 ) {
        # Evt. tekst paa punkter saettes lidt nede til hoejre
        text(x,y,txt,adj=c(-.75,.75),cex=fex)
      }
   }  # if ( add == FALSE )



   if ( RTS == "add" )   {
      # Lav alle mulige additive kombinatinoner af data
      # Find foerst randen i en fdh
      idx <- sort(x, index.return=TRUE)
      effektive <- rep(NA, length(x))
      j <- 1
      prev <- idx$ix[1]
      effektive[j] <- prev
      for ( i in idx$ix )  {
         if ( y[i] > y[prev] )  {
         	j <- j+1
            effektive[j] <- i
            prev <- i
         }
      }
      # Frontier/rand og antal firms i randen
      rand <- effektive[!is.na(effektive)]
      nr <- length(rand)
      # Hvor mange gange en firm optraeder foer vi er uden for plotrammen
      if ( missing(xlim) )  xlim=c(0,1.2*max(x+.001))
      nx <- round(xlim[2]/x[rand])

      # lav aggregerings matrix til alle kombinationer af data
      # Hoejst 5 gentagelser hvis nx er uendelig fordi x[rand] er 0
      if ( is.infinite(nx) )  nx <- 5
      M <- matrix(NA, nrow=prod(1+nx), ncol=nr)
      M[,1] <- rep(0:nx[1], prod(1+nx[-1]))
      if ( nr>1) for ( j in 2:nr )  {
         M[,j] <- as.integer(gl(1+nx[j], prod(1+nx[1:(j-1)]))) -1
      }
      # Drop foerste raekke med bar nuller
      M <- M[-1,]
      # Lav saa alle kombinationerne
      x <- M %*% as.matrix(x[rand])
      y <- M %*% as.matrix(y[rand])
      # Herefter kan plottet laves som var det fdh, bemaerk at firms
      # er allerede afsat som punkter saa kombinationerne kommer ikke
      # til at optraede som punkter
      RTS <- "fdh"
   }



   if ( RTS == "fdh+" )  {
      dea.plot.fdhPlus(x,y,param, ...)
   }  # "fdh+"




   if ( ORIENTATION == "in" ) {
      # inputkravmaengde, input afstandsfunktion
      hpts=chull(c(x,min(x),2*max(x)),c(y,2*max(y),min(y)))
      if ( RTS != "fdh" ) {
         lines(x[hpts],y[hpts],...)
      }
      if ( RTS == "fdh" )  {
         # tegn linjer for fdh
         idx <- sort(x,index.return=T)
         prev <- idx$ix[1];
         for ( i in idx$ix )  {
            if ( y[i] < y[prev] ) {
               lines(x[c(prev,i,i)],y[c(prev,prev,i)],...)
               prev <- i
            }
         }
      }
      x1fmax <- max(x[hpts],na.rm=T)
      x2fmax <- max(y[hpts],na.rm=T)
      lines(c(x1fmax,2*max(x)), c(min(y), min(y)),...)
      lines(c(min(x), min(x)), c(x2fmax,2*max(y)),...)
   } else if (ORIENTATION == "out") {
      # produktionsmulighedsomraade for output, output afstandsfunktion
      hpts=chull( c(x,0,0,max(x)) , c(y,0,max(y),0) )
      # For at vaere sikre paa at alle linjer tegnes laves en lukket graf
      # eller kan der opstaa et hul i fronten
      hpts <- c(hpts, hpts[1])
      if ( RTS != "fdh" ) {
         lines(x[hpts],y[hpts],...)
         # Problem hvis mindste x er 0 ved max(y) for saa bliver
         # min(x(hpts)) ikke 0 da det omtalte punkt ikke er i hpts som
         # del af x, men det ekstra punkt (0,max(y)).
         lines(c(0, min(x[hpts],na.rm=T) ), 
               c(max(y),max(y[hpts],na.rm=T) ),...)
         lines( c(max(x),max(x)), c(0,min(y[hpts],na.rm=T)),...)
      }
      if ( RTS == "fdh" )  {
         # tegn linjer for fdh
         idy <- sort(y,index.return=T,decreasing=T)
         lines(c(0,x[idy$ix[1]]),c(idy$x[1],idy$x[1]),...)
         prev <- match(max(x),x)
         lines(c(max(x),max(x)),c(0,y[prev]),...)
         prev <- idy$ix[1];
         for ( i in idy$ix[-1] )  {
            if ( x[i] > x[prev] ) {
               lines(x[c(prev,prev,i)],y[c(prev,i,i)],...)
               prev <- i
            }
         }
      }
   } else {  # ORIENTATION == "in-out"
      # Et input og et output, "normal produktionsfunktion"
      # Foerst findes det konvekse hul af punkterne og linjer 
      # for yderpunkter tegnes.
      # Punkterne udvides med noget stoerre end max for at tegne linjer
      # der peger mod uendelig (free disposability)
      if ( RTS == "crs" ) {
         # crs
         abline(0, max(y/x),...)
      } else if (RTS=="vrs" | RTS=="fdh")  {
         # vrs
         # hpts=chull(c(x,min(x),max(x)+2),c(y,0,max(y)))
         hpts=chull( c(x,min(x),2*max(x),2*max(x)),c(pmax(y,0),0,max(y),0) )
         lines(c(min(x),min(x)), c(0,min(y[hpts],na.rm=T)),...)
         lines(c(max(x[hpts],na.rm=T),2*max(x)),c(max(y),max(y)),...)
      } else if (RTS=="drs") {
         # vrs og (0,0), dvs. drs
         hpts=chull(c(x,0,2*max(x),2*max(x)),c(pmax(y,0),0,0,max(y)))
         lines(c(0,min(x[hpts],na.rm=T)), c(0,min(y[hpts],na.rm=T)),...)
         lines(c(max(x[hpts],na.rm=T),2*max(x)),c(max(y),max(y)),...)
      } else if (RTS=="irs") {
         # vrs plus infinity
         hpts=chull( c(x ,min(x),2*max(x),          2*max(x)),
                     c(pmax(y,0),0,     2*max(x)*max(y/x),0) )
         # get the unit with the largest slope
         id <- which.max(y/x)
         lines(c(min(x),min(x)), c(0,min(y[hpts],na.rm=T)),...)
         lines(c(x[id],  2*max(x)),
               c(y[id], 2*max(x)*y[id]/x[id]),...)
      }
      # For vrs, drs and irs draw the lines between the 
      # extreme points in the convex hull
      if ( RTS == "vrs" | RTS == "drs" | RTS == "irs" ) 
           lines(x[hpts],y[hpts],...)
      if ( RTS == "fdh" | RTS == "add" )  {
         # tegn linjer for fdh
         idx <- sort(x,index.return=TRUE)
         prev <- idx$ix[1];
         for ( i in idx$ix )  {
            if ( y[i] > y[prev] ) {
               lines(x[c(prev,i,i)],y[c(prev,prev,i)],...)
               prev <- i
            }
         }
      }
      if ( RTS == "irs2" )  {
         # Lines for increasing returns to scale, irs2
         # Plot first vertical part
         lines(c(min(x),min(x)), c(0,max(y[x==min(x)],na.rm=T)),...)

         # Plot the line for input going to infinity        
         slopes <- y/x
         idmax <- which.max(y/x)
         lines(c(x[idmax],20*x[idmax]),
               c(y[idmax],20*slopes[idmax]*x[idmax]),...)
         
         # find inputs in increasing order that have increasing slopes
         sx <- sort(x,index.return=T)
         id <- which(x %in% min(x),y)
         idy <- id[which.max(y[id])]
         slope <- y[idy]/x[idy]
         hpts <- idy
         for( i in 1:length(x) )  {
            if ( slopes[sx$ix[i]] > slope ) {
               # A higher slope found, change level in technology
               hpts <- c(hpts,sx$ix[i])
               slope <- slopes[sx$ix[i]]
            }
         }
         # print(hpts)

         # Plot the lines for increasing slopes, remember the jumps
         if ( length(hpts) > 1 )  {
         for ( i in 1:(length(hpts)-1) )  {
            lines(c(x[hpts[i]],x[hpts[i+1]],x[hpts[i+1]]),
             c(y[hpts[i]],slopes[hpts[i]]*x[hpts[i+1]],y[hpts[i+1]]),...)
         } }
      } # RTS=="irs2"
   }
}



dea.plot.frontier <- function(x, y, RTS="vrs",...)
{
        dea.plot(x, y, RTS=RTS, ORIENTATION="in-out",...)
}

dea.plot.isoquant <- function(x1, x2, RTS="vrs",...)
{
        dea.plot(x1, x2, RTS=RTS, ORIENTATION="in",...)
}

dea.plot.transform <- function(y1, y2, RTS="vrs",...)
{
        dea.plot(y1, y2, RTS=RTS, ORIENTATION="out",...)
}

