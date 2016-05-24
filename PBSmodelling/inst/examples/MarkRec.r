# R Code to run the mark-recovery model

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

# ************************************************************
# Functions required for the GUI
# ************************************************************

# Initialize and compile the WinBUGS model

modCompile <- function() {
  getWinVal(scope="L");              # get data from the GUI
  mdat <- list(M=M,S=S,R=R,eps=eps); # collect input data
  bugsData(mdat,"MarkRecDat.txt");   # write data in bugs format
  modelCheck("MarkRecMod.txt");      # check model syntax
  modelData("MarkRecDat.txt");       # load current data
  modelCompile(nc);                  # compile with nc chains
  modelGenInits();                   # generate randoms inits
  samplesSet(c("p","N"));            # set parameters to monitor
  Nest <- M*S/R;
  setWinVal(list(Nest=Nest));        # update estimate in GUI
}

# Update the model and save complete history in global "MRhist"

modUpdate <- function() {
  getWinVal(scope="L");
  modelUpdate(clen,cthin);
  MRhist <- as.data.frame( samplesHistory("*",beg=0,plot=FALSE) ); tput(MRhist)
  ctot <- dim(MRhist)[1];    # total length so far
  setWinVal(list(ctot=ctot)); par(ask=FALSE); }
  
# Functions to report results

modHist <- function() {
  getWinVal(scope="L"); resetGraph();
  samplesHistory("*",beg=s1-1,end=s2-1,thin=sthin,
    mfrow=c(2,1),ask=FALSE); };

modDens <- function() {
  getWinVal(scope="L"); resetGraph();
  samplesDensity("*",beg=s1-1,end=s2-1,thin=sthin,
    mfrow=c(2,1),ask=FALSE); };

modACF <- function(chn=1) { # default chain 1
  getWinVal(scope="L"); resetGraph();
  samplesAutoC("*",chn,beg=s1-1,end=s2-1,thin=sthin,
    mfrow=c(2,1),ask=FALSE); };

modPairs <- function() {
  getWinVal(scope="L"); resetGraph();
  i1 <- max(s1,1); i2 <- min(s2,ctot); # ensure valid range
  idx <- seq(i1,i2,by=sthin);
  tget(MRhist)
  par(ask=FALSE); pairs(MRhist[idx,],pch=19,cex=0.5); };

# ************************************************************
# Load libraries and start the GUI
# ************************************************************

if (!require(BRugs, quietly=TRUE)) stop("The `BRugs` package is required for this example")
if (!require(PBSmodelling, quietly=TRUE)) stop("The `PBSmodelling` package is required for this example")
createWin("MarkRecWin.txt")

}) # end local scope
