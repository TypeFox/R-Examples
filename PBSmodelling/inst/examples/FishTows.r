# R Code to run the tow simulation model

# Demo for PBS Modelling, dependent on PBS Mapping
# Simulation to generate a random set of tows
# within a given square region. See FishTowsDoc.txt.

# ***********************************
# Functions for dealing with tow data
# ***********************************

# By definition, a set of tows (called a TowSet) is a PolySet of polylines with
# only two points, the start and end of each tow.  A TowSet contains only PIDs,
# not SIDs, and has two POS values 1 and 2 for each tow.


# Function to generate a random set of tows within a given square region.
# Input:  n = number of tows
#         s = side length of a square from the origin
# Output: a TowSet

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

makeTows <- function(n=10,s=100) {
   X <- runif(2*n,0,s); Y <- runif(2*n,0,s);
   ieven <- 2*(1:n); iodd <- ieven-1;
   PID <- rep(0,2*n); POS <- PID;
   PID[iodd] <- 1:n; PID[ieven] <- 1:n;
   POS[iodd] <- 1; POS[ieven] <- 2;
   tset <- data.frame(PID,POS,X,Y)
   attr(tset,"projection") <- 1;
   attr(tset,"class") <- c("PolySet",attr(tset,"class"));
   invisible(tset); };

# Function to expand tows in a TowSet to rectangles with a specified width.
# Warning: no check that the TowSet is valid
# Input:  tset = a specified TowSet
#         w = width of a tow
# Output: a PolySet of rectangles

expandTows <- function(tset,w) {
   tpos <- tset$POS; tproj <- attr(tset,"projection");
   x1 <- tset$X[tpos==1]; x2 <- tset$X[tpos==2];
   y1 <- tset$Y[tpos==1]; y2 <- tset$Y[tpos==2];
   r <- sqrt( (x2-x1)^2 + (y2-y1)^2 ); # tow lengths
   xd <- w*(y2-y1)/(2*r);
   yd <- - w*(x2-x1)/(2*r);
   x1a <- x1 + xd; x1b <- x1 - xd;
   x2a <- x2 - xd; x2b <- x2 + xd;
   y1a <- y1 + yd; y1b <- y1 - yd;
   y2a <- y2 - yd; y2b <- y2 + yd;
   tdat <- data.frame(EID=1:length(x1),
      x1a, x1b, x2a, x2b, y1a, y1b, y2a, y2b);
   attr(tdat,"projection") <- tproj;
   towP <- convDP(tdat, xColumns=c("x1a","x1b","x2a","x2b"),
      yColumns=c("y1a","y1b","y2a","y2b"));
   # attr(towP,"projection") <- 1;
   invisible(towP); };

# *************************
# Extensions to PBS Mapping
# *************************

# Find holes in a PolySet (based on summary.PolySet)

findHoles <- function(PSet) {
   if (!is.element("SID", names(PSet))) return(NULL);
   IDX <- .createIDs(PSet, cols = c("PID", "SID")); # numeric IDs
   idxFirst <- which(!duplicated(IDX));
   idxLast <- c((idxFirst - 1)[-1], length(IDX));
   hole <- (PSet$POS[idxFirst] > PSet$POS[idxLast]);
   data.frame(
      PID=PSet$PID[idxFirst],
      SID=PSet$SID[idxFirst],
      row=idxFirst,
      IDX=IDX[idxFirst], hole=hole); };

# Divide a PolySet so that each PID contains a single main polygon
# with holes in corresponding SIDs

divPolys <- function(PSet) {
   if (!is.element("SID", names(PSet))) return(PSet);
   IDX <- .createIDs(PSet, cols = c("PID", "SID")); # numeric IDs
   idxFirst <- which(!duplicated(IDX));
   nrec   <- length(IDX);               # number of records
   npoly  <- length(idxFirst);          # number of polygons (main & hole)
   idxLast<- c((idxFirst - 1)[-1], nrec);
   hole   <- (PSet$POS[idxFirst] > PSet$POS[idxLast]);
   plen   <- idxLast - idxFirst + 1;    # length (records) in each polygon
   main   <- !hole;
   nmain  <- sum(main);                 # number of main polygons
   idxm   <- idxFirst[main];            # starting indices for mains
   mainL  <- c(idxm[-1],nrec+1) - idxm; # lengths of mains
   newPID <- rep(1:nmain,mainL);
   newSID <- rep(1:npoly,plen);
   for (i in (1:nmain)) {					 # number SIDs from 1 for each main
      irow <- idxm[i];                  # initial row index (main i)
      prng <- 1:mainL[i];
      irng <- irow + prng - 1;          # range of indices (main i)
      p1   <- newSID[irow];
      newSID[irng] <- newSID[irng] - p1 + 1; };
   POut <- PSet; POut$PID <- newPID; POut$SID <- newSID;
   POut; };

genTows <- function() {
   getWinVal(scope="L");
   n <- nTow; w <- wTow; s <- sTow;
   towLines <- makeTows(n,s); tput(towLines)
   towPolys <- expandTows(towLines,w); tput(towPolys)
   towAll   <- joinPolys(towPolys, operation="UNION")
   towAll   <- divPolys(towAll); tput(towAll)
   towSumm <- summary(towAll);
   nCont   <- towSumm$contours; nHole <- towSumm$holes;
   nPoly   <- nCont - nHole;
   nVert   <- towSumm$records;
   swpArea <- sum(calcArea(towPolys,rollup=1)$area);
   impArea <- sum(calcArea(towAll,rollup=1)$area);
   totArea <- s^2;
   meanL   <- swpArea/(n*w);
   setWinVal(c(meanL=signif(meanL,6),nPoly=nPoly,nHole=nHole,nVert=nVert,
      swpArea=signif(swpArea,8),impArea=signif(impArea,8),totArea=totArea));
   invisible(); };

plotTow <- function(act=NULL) {
   resetGraph(); getWinVal(scope="L");
   s <- sTow; eps <- s/50;
   # range slightly larger than (0,s)
   rng0 <- c(-eps,s+eps);
   # larger range for rectangles
	tget(towPolys)
   rng1 <- range(c(-eps,towPolys$X,towPolys$Y,s+eps));
   if (is.null(act)) act <- getWinAct()[1];
   side <- 6 # side of 1 frame in inches for EMF
   if (wmf) win.metafile(filename=paste("FishTows",act,".emf",sep=""),pointsize=12,
       width=ifelse(act=="C" && cmode=="L",2*side,side), height=ifelse(act=="C" && cmode=="P",2*side,side));
   par(mar=c(2,2,0.1,0.5),mgp=c(2,.25,0),las=1);
   if (act=="L") {
   	tget(towLines)
      plotMap(towLines,type="n",xlim=rng0,ylim=rng0,plt=NULL,cex=1.5,xlab="",ylab="");
     addLines(towLines); };
   if (act=="T")
      plotMap(towPolys,col="red",xlim=rng1,ylim=rng1,plt=NULL,cex=1.5,xlab="",ylab="");
   if ((act=="U") | (act=="C")) {
      myCol <- c("red","green","blue","yellow","orange","violet","cyan",
         "firebrick","deeppink","brown","chocolate","coral","azure");
      myCol <- rep(myCol,100); # up to 1300 colors
      tget(towAll)
      towID <- unique(towAll$PID); nID <- length(towID);
      towCol <- myCol[1:nID];
      towDF <- data.frame(PID=towID,col=towCol);
      towDF[[2]] <- as.character(towDF[[2]]);
      towProps <- as.PolyData(towDF);
      if (act=="U")
         plotMap(towAll,polyProps=towProps,xlim=rng1,ylim=rng1,colHoles="white",plt=NULL,cex=1.5,xlab="",ylab="");
      if (act=="C") {
         if (cmode=="P")  par(mfrow=c(2,1),mar=c(1.2,1.4,0.1,0.5),mgp=c(2,.1,0),las=1);
         if (cmode=="L")  par(mfrow=c(1,2),mar=c(1.2,1.4,0.1,0.5),mgp=c(2,.1,0),las=1);
         plotMap(towPolys,col=myCol[1],xlim=rng1,ylim=rng1,
            xlab="",ylab="",plt=NULL,cex=1.5);
         plotMap(towAll,polyProps=towProps,xlim=rng1,ylim=rng1,
            xlab="",ylab="",plt=NULL,cex=1.5,colHoles="white");
   } };
   if (wmf) dev.off();
   invisible(); };

# ******************
# Initialize the GUI
# ******************

if (!require(PBSmapping, quietly=TRUE)) stop("The `PBSmapping` package is required for this example")
if (!require(PBSmodelling, quietly=TRUE)) stop("The `PBSmodelling` package is required for this example")
createWin("FishTowsWin.txt")

}) # end local scope
