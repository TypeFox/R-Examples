# Figures for PBSmapping examples (last modified: 2015-02-17)
#------------------------------------------------------------
# Historical values for compatibilityy with S-Plus (defunct)
.PBSdot <- 3; .PBSdash <- 2
.PBSclr <- function(){
   PBSclr = list(             black=c(0,0,0),
      sea=c(224,253,254),     land=c(255,255,195),      red=c(255,0,0),
      green=c(0,255,0),       blue=c(0,0,255),          yellow=c(255,255,0),
      cyan=c(0,255,255),      magenta=c(255,0,255),     purple=c(150,0,150),
      lettuce=c(205,241,203), moss=c(132,221,124),      irish=c(54,182,48),
      forest=c(29,98,27),     white=c(255,255,255),     fog=c(223,223,223) )
      PBSclr <- lapply(PBSclr,function(v) {rgb(v[1],v[2],v[3],maxColorValue=255) })
      return(PBSclr) }

.PBSfig01 <- function() { #  World UTM Zones
   clr <- .PBSclr()
   data(worldLL,nepacLL,envir=sys.frame(sys.nframe()))
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(worldLL, ylim=c(-90, 90), bg=clr$sea, col=clr$land, tck=-0.023,
           mgp=c(1.9, 0.7, 0), cex=1.2, plt=c(.08,.98,.08,.98))
   # add UTM zone boundaries
   abline(v=seq(-18, 360, by=6), lty=1, col=clr$red)
   # add prime meridian
   abline(v=0, lty=1, lwd=2, col=clr$black)
   # calculate the limits of the 'nepacLL' PolySet
   xlim <- range(nepacLL$X) + 360
   ylim <- range(nepacLL$Y)
   # create and then add the 'nepacLL' rectangle
   region <- data.frame(PID=rep(1,4), POS=1:4, X=c(xlim[1],xlim[2],xlim[2],xlim[1]),
                        Y=c(ylim[1],ylim[1],ylim[2],ylim[2]))
   region <- as.PolySet(region, projection="LL")
   addPolys(region, lwd=2, border=clr$blue, density=0)
   # add labels for some UTM zones
   text(x=seq(183.2, by=6, length=9), y=rep(85,9), adj=0.5, cex=0.65, label=1:9)
   box() }

.PBSfig02 <- function() {  # nepacLL UTM Zones in LL Space
   clr <- .PBSclr(); dot <- .PBSdot
   data(nepacLL,envir=sys.frame(sys.nframe()))
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(nepacLL, col=clr$land, bg=clr$sea, tck=-0.014,
           mgp=c(1.9,0.7,0), cex=1.2, plt=c(.08,.98,.08,.98))
   # add lines separating UTM zones
   utms <- seq(-186, -110, 6)
   abline(v=utms, col=clr$red)
   # add the central meridian of zone 6
   abline(v=-147, lty=dot, col=clr$black)
   # create and then add labels for the UTM zones
   cutm <- diff(utms) / 2
   nzon <- length(cutm)
   cutm <- cutm + utms[1:nzon]
   text(cutm,rep(50.75,nzon),c(60,1:(nzon-1)),cex=1.3,col=clr$red)
   box()  }

.PBSfig03 <- function() {  # nepacLL UTM Zones in UTM Space
   clr <- .PBSclr(); dot <- .PBSdot
   data(nepacLL,envir=sys.frame(sys.nframe()))
   zone  <- 6;  xlim  <- range(nepacLL$X);  ylim <- range(nepacLL$Y)
   utms  <- seq(-186,-110,6)  #'utms' vector for creating PolySet and EventData below
   # create UTM zones
   lutms <- data.frame(PID=rep(1:length(utms), each=2),
               POS=rep(c(1,2), times=length(utms)), X=rep(utms,each=2),
               Y = rep(c(ylim[1], ylim[2]), times=length(utms)))
   lutms <- as.PolySet(lutms, projection="LL", zone=zone)
   lutms <- thickenPolys(lutms, tol=25, close=FALSE)
   uutms <- convUL(lutms)
   # create label locations (central meridians)
   lcms  <- data.frame(EID=1:(length(diff(utms)/2)),
               X=utms[1:(length(utms)-1)]+diff(utms)/2,
               Y=rep(50.75, length(diff(utms)/2)))
   lcms  <- as.EventData(lcms, projection="LL", zone=zone)
   ucms  <- convUL(lcms)
   nepacUTM <- nepacLL; attr(nepacUTM,"zone") <- zone  # convert to correct zone
   nepacUTM <- convUL(nepacUTM)
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(nepacUTM, col=clr$land, bg=clr$sea, tck=-0.017,
           mgp=c(1.9,0.7,0), cex=1.0, plt=c(0.07,0.97,0.07,0.98))
   addLines(uutms, col=clr$red)
   lines(x=c(500, 500),y=c(4100,7940),lty=dot,col=clr$black)
   text(ucms$X,ucms$Y,c(60,1:(length(utms)-2)),cex=1.3,col=clr$red)
   box()  }

.PBSfig04 <- function() {  # thinPolys on Vancouver Island
   clr <- .PBSclr();
   data(nepacLL,envir=sys.frame(sys.nframe()))
   par(mfrow=c(1,2),omi=c(0,0,0,0)) #------Plot-the-figure------
   vi     <- nepacLL[nepacLL$PID==33,]
   xlim   <- range(vi$X) + c(-0.25, 0.25); ylim <- range(vi$Y) + c(-0.25, 0.25)
   # plot left figure (normal Vancouver Island)
   plotMap(vi, xlim, ylim, col=clr$land, bg=clr$sea, tck=-0.028,
           mgp=c(1.9,0.7,0), cex=1.0, plt=c(0.14,1.00,0.07,0.97))
   text(x=xlim[2]-0.5, y=ylim[2]-0.3, "A", cex=1.6)
   # plot right figure (thinned Vancouver Island)
   plotMap(thinPolys(vi, tol=10), xlim, ylim, col=clr$land, bg=clr$sea,
           tck=c(-0.028, 0), tckLab=c(TRUE, FALSE),
           mgp=c(1.9, 0.7, 0), cex=1.0, plt=c(0.00, 0.86, 0.07, 0.97))
   text(x=xlim[2]-0.5, y=ylim[2]-0.3, "B", cex=1.6)
   box() }

.PBSfig05 <- function() {  # joinPolys on Crescents
   clr <- .PBSclr(); dash <- .PBSdash
   radius <- c(5, 4)                 # two radii of the circles
   size   <- abs(diff(radius)) + 0.1 # size of crescent
   shiftB <- 3.5                     # distance to shift second crescent
   pts    <- 120                     # points in outer circle
   cex    <- 1.0                     # character expansion for labels
   off    <- 1.2                     # panel label offset
   xlim   <- c(0, radius[1]*2 + shiftB) + c(-1,1)
   ylim   <- c(0, radius[1]*2) + c(-2,1)
   Mmin   <- .10 # minimum OMI
   Rdin   <- par()$din[2]/par()$din[1]
   Rfig   <- (3*diff(ylim))/(2*diff(xlim))
   if (Rdin > Rfig) {
      width  <- par()$din[1] - 2 * Mmin
      height <- width * (3*diff(ylim))/(2*diff(xlim))
      Mmax   <- (par()$din[2] - height) / 2
      parOmi <- c(Mmax,Mmin,Mmax,Mmin) }
   else {
      height <- par()$din[2] - 2 * Mmin
      width  <- height * (2*diff(xlim))/(3*diff(ylim))
      Mmax   <- (par()$din[1] - width) / 2
      parOmi <- c(Mmin,Mmax,Mmin,Mmax) }
   polyA  <- list()
   for (i in 1:length(radius)) {
      polyA[[i]] <- as.PolySet(data.frame(PID=rep(1,pts), POS = 1:pts,
         X =radius[i]*cos(seq(0, 2*pi, len=pts)),
         Y =radius[i]*sin(seq(0, 2*pi, len=pts))), projection = 1)
      polyA[[i]][, c("X","Y")] <- polyA[[i]][, c("X","Y")] + radius[i] }
   # centre B within A
   polyA[[2]][,c("X","Y")] <- polyA[[2]][,c("X","Y")] + (radius[1]-radius[2])
   # shift B right
   polyA[[2]]$X <- polyA[[2]]$X + size
   # create 'polysA' and 'polysB'
   polyA  <- as.PolySet(joinPolys(polyA[[1]], polyA[[2]], operation="DIFF"), projection=1)
   polyB  <- polyA
   polyB$X<- abs(polyB$X - (radius[1] * 2)) + shiftB
   par(mfrow=c(3,2),mai=c(0,0,0,0),omi=parOmi) #------Plot-the-figure------
   lab    <- list()
   lab$text <- c("Polygon A", "Polygon B", "A \"INT\" B","A \"UNION\" B",
                 "A \"DIFF\" B", "A \"XOR\" B")
   lab$cex <- rep(cex, 6);  lab$x <- rep(mean(xlim), 6);  lab$y <- rep(-0.8, 6)
   # panel A: polyA
   plotMap(polyA,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,col=clr$red,plt=NULL)
   text(lab$text[1], x=lab$x[1], y=lab$y[1], cex=lab$cex[1])
   text(xlim[1]+off, ylim[2]-off, "A", cex=1.6);  box()
   # panel B: polyB
   plotMap(polyB,xlim=xlim,ylim=ylim,xlab="",ylab="",axes=FALSE,col=clr$blue,plt=NULL)
   text(lab$text[2], x=lab$x[2], y=lab$y[2], cex=lab$cex[2])
   text(xlim[1]+off, ylim[2]-off, "B", cex=1.6);  box()
   # panels C to F
   ops    <- c(NA, NA, "INT", "UNION", "DIFF", "XOR")
   cols   <- c(NA, NA, clr$red, clr$purple, clr$red, clr$magenta)
   panel  <- c(NA, NA, "C", "D", "E", "F")
   for (i in 3:6) {
      plotMap(NULL,xlim=xlim,ylim=ylim,projection=1,xlab="",ylab="",axes=FALSE,plt=NULL)
      addPolys(polyA, border=clr$red, lty=dash)
      addPolys(polyB, border=clr$blue, lty=dash)
      addPolys(joinPolys(polyA, polyB, operation=ops[i]), col=cols[i])
      text(lab$text[i], x=lab$x[i], y=lab$y[i], cex=lab$cex[i])
      text(xlim[1]+off, ylim[2]-off, panel[i], cex=1.6);  box();  } }

.PBSfig06 <- function() {  # contourLines in Queen Charlotte Sound
   clr <- .PBSclr(); 
   data(nepacLL,bcBathymetry,envir=sys.frame(sys.nframe()));
   isob   <- contourLines(bcBathymetry, levels=c(250, 1000))
   p      <- convCP(isob)
   attr(p$PolySet,"projection") <- "LL"
   p$PolyData$col <- rep(c(clr$red, clr$green, clr$blue, clr$yellow,
      clr$cyan, clr$magenta, clr$fog), length=nrow(p$PolyData))
   xlim   <- c(-131.8382, -128.2188)
   ylim   <- c(50.42407, 53.232476)
   region <- clipPolys(nepacLL, xlim=xlim, ylim=ylim)
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #-----Plot-the-figure------
   plotMap(region, xlim=xlim, ylim=ylim, col=clr$land, bg=clr$sea, tck=-0.02,
        mgp=c(2,.75,0), cex=1.2, plt=c(.08,.98,.08,.98))
   addLines(p$PolySet, polyProps=p$PolyData, lwd=3)
   box()  }

.PBSfig07 <- function() {  # towTracks from Longspine Thornyhead Survey
   clr <- .PBSclr();
   data(nepacLL,towTracks,towData,envir=sys.frame(sys.nframe()));
   # add a colour column 'col' to 'towData'
   pdata  <- towData;  pdata$Z <- pdata$dep
   pdata  <- makeProps(pdata, breaks=c(500,800,1200,1600), "col",
                       c(clr$black, clr$red, clr$blue))
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(nepacLL, col=clr$land, bg=clr$sea, xlim=c(-127.8,-125.5), ylim=c(48,49.8),
      tck=-0.01, mgp=c(2,.5,0), cex=1.2, plt=c(.08,1,.08,.98))
   addLines(towTracks, polyProps=pdata, lwd=3)
   # right-justify the legend labels
   temp <- legend(x=-127.6, y=48.4, legend=c(" "," "," "), lwd=3, bty="n",
      text.width=strwidth("1200-1600 m"), col=c(clr$black,clr$red,clr$blue))
   text(temp$rect$left+temp$rect$w, temp$text$y,
      c("500-800 m", "800-1200 m", "1200-1600 m"), pos=2)
   text(temp$rect$left+temp$rect$w/2,temp$rect$top,pos=3,"LTS Survey Tracks");
   text(-125.6,49.7,"Vancouver\nIsland",cex=1.2,adj=1)
   box()  }

.PBSfig08 <- function() {  # calcArea of the Southern Gulf Islands
   clr <- .PBSclr(); 
   data (nepacLLhigh,envir=sys.frame(sys.nframe()))
   xlim   <- c(-123.6, -122.95); ylim <- c(48.4, 49); zone <- 9
   # assign 'nepacLLhigh' to 'nepacUTMhigh' (S62) and change to UTM coordinates
   nepacUTMhigh <- nepacLLhigh;  attr(nepacUTMhigh,"zone" ) <- zone
   nepacUTMhigh  <- convUL(nepacUTMhigh)
   # convert limits to UTM
   temp   <- data.frame(PID=1:4,POS=rep(1,4),X=c(xlim,xlim),Y=c(ylim,rev(ylim)))
   temp   <- convUL(as.PolySet(temp, projection="LL", zone=zone))
   xlim   <- range(temp$X); ylim <- range(temp$Y)
   # prepare areas
   isles  <- clipPolys(nepacUTMhigh,xlim,ylim)
   areas  <- calcArea(isles);
   # PIDs and labels for Gulf Islands
   bigPID <- areas[rev(order(areas$area)),][c(2:4,6:8),"PID"];
   labelData <- data.frame(PID = bigPID, 
      label=c("Saltspring","San Juan","Galiano","Saturna","N Pender","Mayne"))
   labelData <- merge(labelData, areas, all.x=TRUE)
   labelData$label <- paste(as.character(labelData$label),
      round(labelData$area), sep="\n")
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(isles, col=clr$land, bg=clr$sea, tck=-.010,
      mgp=c(1.9,.7,0), cex=1, plt=c(.07,.98,.07,.98))
   # add the highlighted Gulf Islands
   bigisles <- isles[is.element(isles$PID,labelData$PID),]
   addPolys(bigisles,col=clr$yellow)
   labXY  <- calcCentroid(isles)
   labXY$Y<- labXY$Y + 2               # centre vertically
   labelData <- merge(labelData, labXY, all.x = TRUE)
   attr(labelData,"projection") <- "UTM"
   addLabels(labelData, placement="DATA", cex=1.25)
   text(898,5385,"Vancouver Island",adj=0, cex=1.25)
   text(925,5435,"Strait of Georgia",adj=0, cex=1.25)  }

.PBSfig09 <- function() {  # combineEvents in Queen Charlotte Sound
   clr <- .PBSclr(); 
   data(nepacLL,surveyData,envir=sys.frame(sys.nframe()));
   events <- surveyData
   xl     <- c(-131.8, -127.2);  yl <- c(50.5, 52.7)
   # prepare EventData; clip it, omit NA entries, and calculate CPUE
   events <- events[events$X >= xl[1] & events$X <= xl[2] &
                    events$Y >= yl[1] & events$Y <= yl[2], ]
   events <- na.omit(events)
   events$cpue <- events$catch/(events$effort/60)
   # make a grid for the Queen Charlotte Sound
   grid   <- makeGrid(x=seq(-131.6,-127.6,.1), y=seq(50.6,52.6,.1),
                    projection="LL", zone=9)
   # locate EventData in grid
   locData<- findCells(events, grid)
   events$Z <- events$cpue
   pdata  <- combineEvents(events, locData, FUN=mean)
   brks   <- c(0,50,300,750,1500,25000); lbrks <- length(brks)
   cols   <- c(clr$lettuce, clr$moss, clr$irish, clr$forest, clr$black)
   pdata  <- makeProps(pdata, brks, "col", cols)
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(nepacLL, col=clr$land, bg=clr$sea, xlim=xl, ylim=yl, tck=-0.015,
           mgp=c(2,.5,0), cex=1.2, plt=c(.08,.98,.08,.98))
   addPolys(grid, polyProps=pdata)
   for (i in 1:nrow(events)) {
      # plot one point at a time for clarity
      points(events$X[i], events$Y[i], pch=16,cex=0.50,col=clr$white)
      points(events$X[i], events$Y[i], pch=1, cex=0.55,col=clr$black) }
   yrtxt  <- paste("(",min(events$year),"-",
                substring(max(events$year),3),")",sep="")
   text(xl[1]+.5,yl[2]-.1,paste("POP Surveys",yrtxt),cex=1.2,adj=0)
   # add a legend; right-justify the legend labels
   temp <- legend(x=xl[1]+.3, y=yl[1]+.7, legend = rep(" ", 5),
              text.width=strwidth("1500 - 25000"), bty="n", fill=cols)
   text(temp$rect$left + temp$rect$w, temp$text$y, pos=2,
        paste(brks[1:(lbrks-1)],brks[2:lbrks], sep=" - "))
   text(temp$rect$left+temp$rect$w/2,temp$rect$top,pos=3,"CPUE (kg/h)",cex=1);  }

.PBSfig10 <- function() {  # Pythagoras' Theorem Visualized
   clr <- .PBSclr(); 
   data(pythagoras,envir=sys.frame(sys.nframe()))
   # create properties for colouring the polygons
   pythProps <- data.frame(PID=c(1, 6:13, 4, 15, 3, 5, 2, 14),
                   Z=c(rep(1, 9), rep(2, 2), rep(3, 2), rep(4, 2)))
   pythProps <- makeProps(pythProps, c(0, 1.1, 2.1, 3.1, 4.1), "col",
                   c(clr$blue, clr$red, clr$yellow, clr$green))
   par(mfrow=c(1,1),omi=c(0,0,0,0)) #------Plot-the-figure------
   plotMap(pythagoras, plt=c(.01,.99,.01,.95), lwd=2,
      xlim=c(.09,1.91), ylim=c(0.19,2.86), polyProps=pythProps,
      axes=FALSE, xlab="", ylab="",
      main=bquote(paste("Pythagoras' Theorem: ", a^2 +  b^2 == c^2)))
   text(x = 0.1, y = 1.19, adj=0, "Proof:")
   text(x = 0.1, y = 1.10, adj=0,
      bquote(paste((a + b)^2 == 4, " triangles ", + a^2 + b^2 == 4,
                   " triangles ", + c^2)))
   labels <- data.frame(X=c(1.02,1.66,0.65),Y=c(1.50,2.20,2.76),label=c("a","b","c"))
   text(labels$X, labels$Y, as.character(labels$label), cex=1.2)
   text(1.03, 1.81, bquote(a^2), cex=1.2, col=clr$black)
   text(1.43, 2.21, bquote(b^2), cex=1.2, col=clr$black)
   text(0.87, 2.46, bquote(c^2), cex=1.2, col=clr$black)  }

.PBSfigs <- function(nfigs=1:10) { # Draw all figures with numbers in nfigs
   #while (!is.null(dev.list())) dev.off(dev.cur())
   for (i in nfigs) {
      figStr <- paste(".PBSfig",ifelse(i<10,"0",""),i,sep="")
      get(figStr)();
      cat(figStr); readline(); }  }
