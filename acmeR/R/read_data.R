#' Construct ACME Sufficient Statistics
#' 
#' Reads in an event-level dataset of carcass placements and searches and
#' constructs a carcass-level and search-level table of sufficient
#' statistics.
#' 
#'@param fname Data, either a string for csv files or a data frame name 
#' @param spec Species subset. Default is empty string.
#' @param blind logial. If TRUE, ensures FT are always unaware
#'  of carcasses
#' @param tz Time Zone. Default is US West Coast
#' @return \code{read.data} returns an invisible list with components:
#' \item{scav}{carcass-level table of removal data}
#' \item{srch}{event-level table or searcher proficiency data}
#' \item{Ik}{summary (count, average, sd) of FT Search Intervals}
#' \item{Sk}{summary (count, average, sd) PFM check intervals}
#' \item{NP.Spec}{number of without a "Placed" event}
#' \item{NP.ID}{number of birds without a "Placed" event}
#' \item{fn}{name of data - parameter \code{fname}}
#' \item{Info}{list of select system information}
#' 
#' @import utils
#' 
#' 
##################################################################
# read.data(): Construct sufficient statistics from input dataset.
#             Read spreadsheet and generate two tables:
#    scav:  Removal data:
#            Id = ID Label
#       Species = Species
#        Placed = Date of placement
#            Lo = Last  date on which carcass is known to be present
#            Hi = First date on which carcass is known to be absent
#    srch: Searcher Proficiency Data:
#            Id = ID Label
#          Date = Dates of searches
#         Found = Discovery result (T=found, F=missed)
#
# Presently ignoring Cal/WEA "birds$AssignedPresence" field and using
# only PFM and FT search results; this changes last-known-date for
# birds: 11  59  62  89 126 127 157 170.  It looks better to me.
#
read.data <- function(fname="acme-sim.csv", 
                       spec="",          
                      blind=TRUE,        # FTs always unaware of carcasses?
                         tz="PST8PDT")   # Time zone (defalt: US West Coast)
{
  
  if(class(fname)=="character"){
    if(!file.exists(ofname<-fname) &&
         !file.exists(fname<-paste(fname,"csv",sep=".")))
      stop("File: \"", ofname, "\" not found");
    input="csv"
    birds <- read.csv(fname, as.is=TRUE, stringsAsFactors=FALSE);
  } else if(class(fname) =="data.frame"){
    input="df"
    birds <- fname
  } else{
   stop("Data fname must be either a file path to a CSV file, a data frame,
        or a matrix.") 
  }
  # Required fields:
  fields    <- c("Date", "ID", "Event", "Found");
  if(any(! (ok <- fields %in% names(birds))))
    stop("Required fields: {", paste(fields[!ok], collapse=", "),
         "} missing from data csv file \"", fname, "\".");
  # Trim down to species requested
  if(!missing(spec) && nchar(spec)) {
    birds  <- birds[ok <- which(birds$Species %in% class2spec(spec)),];
    if(!any(ok))
      stop("File: \"", fname, "\" contains no entries of species: ", spec);
  }
  #
  nEvent <- dim(birds)[1];
  Placed <- grep("^[pP]", birds$Event);   # vector of     Place  indices
  PFM    <- grep("^[cC]", birds$Event);   # vector of PFM Check  indices
  FT     <- grep("^[sS]", birds$Event);   # vector of FT  Search indices
  isPFM  <- 1:nEvent %in% PFM;            # vector of T/F
  isFT   <- 1:nEvent %in% FT;
  isPla  <- 1:nEvent %in% Placed;
  # Optional fields:
  if(! "Species" %in% names(birds)) { birds$Species <- "MISC"; }
  if(! "Time"    %in% names(birds)) {
    birds$Time         <- numeric(length(birds$Event));
    birds$Time[Placed] <-  "8:00:00 AM";  # Offset the times, to
    birds$Time[FT]     <- "12:00:00 PM";  # avoid zero carcass ages
    birds$Time[PFM]    <-  "4:00:00 PM";  # PFM after FT for same-day
  }
  # Cobble together the date+time:
  Days <- birds$Date;
  Time <- birds$Time;
  # ISO 8601 date/time format?
  ISO.d <- grepl("^[0-9]{4,4}-[0-9]{2,2}-[0-9]{2,2} *$", Days);
  ISO.t <- grepl("^ *[0-9]{2,2}:[0-9]{2,2}:[0-9]{2,2} *$", Time);
  USA.d <- grepl("^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4,4} *$", Days);
  USA.t <- grepl("^ *[0-9]{1,2}:[0-9]{2,2}:[0-9]{2,2} *[AaPp][Mm]$", Time);
  
  ISO.days <- all(ISO.d)
  ISO.times <- all(ISO.t)
  USA.days <- all(USA.d)
  USA.times <- all(USA.t)
  if(ISO.days & ISO.times) {
    fmt.dt <- "%Y-%m-%d %H:%M:%S";                # YYYY-mm-dd, 24-hr time
  } else if(ISO.days & USA.times){
    fmt.dt <- "%Y-%m-%d %I:%M:%S %p";             # YYYY-mm-dd, 12-hr time
  } else if(USA.days & USA.times) {
    fmt.dt <- "%m/%d/%Y %I:%M:%S %p";             # mm/dd/YYYY, 12-hr AM/PM
  } else if(USA.days & ISO.times){
    fmt.dt <- "%m/%d/%Y %H:%M:%S";                # mm/dd/YYYY, 24-hr
  } else {
    stop("Expecting date format YYYY-MM-DD (ISO) or MM/DD/YYYY (USA)");
  }
  Date <- strptime(apply(cbind(Days,Time),1,paste,collapse=" "),
                   format=fmt.dt, tz=tz);
  nbirds   <- length(Placed);
  IDs      <- as.character(birds$ID[Placed]);           # List of distinct IDs
  
  #Fix searcher blindness
  if(!blind) {
    for(i in 1:nbirds) {
      srch.i <- isFT   &(birds$ID == IDs[i]);  # FT Search for this bird
      succ.i <- srch.i & birds$Found;          # Found it!
      if(sum(succ.i)>1) {                      # At least twice
        first.i <- min(Date[succ.i]);          # First time; blindness ends.
        other.i <- succ.i & Date > first.i;    # Unblind FT searches
        isFT [other.i] <- FALSE;               # Treat un-blind as TFM check
        isPFM[other.i] <- TRUE;
      }
    }
    PFM <- which(isPFM);
    FT  <- which(isFT);
  }
  
  #Build scavenger dataset
  sID <- sort(IDs);
  if(any(sID[-1]==rev(rev(sID)[-1]))) stop("Duplicate Placement for same ID");
  min.date <- as.POSIXlt("0001-01-01 00:00:00",tz=tz);  # proxy for inf past
  max.date <- as.POSIXlt("9999-12-31 00:00:00",tz=tz);  # proxy for inf future
  scav <- data.frame(stringsAsFactors=FALSE,
                     Species=rep(NA,       nbirds),
                     Id=birds$ID[Placed],
                     Placed=rep(min.date, nbirds),
                     Lo=rep(min.date, nbirds),       # Last  verified presence
                     Hi=rep(max.date, nbirds));      # First verified absence
  nDays <- numeric(nEvent);
  for(i in 1:nbirds) {                 # First, how quickly is it REMoved?
    ok  <- birds$ID==IDs[i];           # Logical vector: this carcass??
    if(any(ok)) {
      scav$Species[i]<- birds$Species[Placed[i]];
      scav$Placed[i] <- Date[Placed[i]];
                                       # Last discovery by PFM or FT:
      for(j in which(ok)) {
         nDays[j]    <- as.numeric(difftime(Date[j],
                                    Date[Placed[i]], units="days"));
      }
      scav$Lo[i]     <- max(Date[ok & birds$Found], Date[Placed[i]]);
                                       # First failure  by PFM (only):
      scav$Hi[i]     <- min(Date[ok & !birds$Found & isPFM &
                                 Date >= scav$Lo[i]],max.date);
    }
  }
  
  #Find carcass age, build interval metrics
  age  <- numeric(nEvent);              # Age in days at event time
  NotPlaced <- numeric(0);              # Anyone never placed?
  for(i in 1:nEvent) {
    arrive <- isPla & (birds$ID == birds$ID[i]);
    if(any(arrive)) {
      age[i] <- as.numeric(difftime(Date[i],Date[arrive],units="days"));
    } else {
      NotPlaced <- c(NotPlaced, i);
      age[i] <- as.numeric(difftime(Date[i],min(Date),units="days"));
    }
  }
  NP.Spec   <- unique(birds$Species[NotPlaced]);
  NP.ID     <- unique(birds$ID[NotPlaced]);
  Cint <- Sint <- numeric(0);          # Crude: Lists of all intervals
  for(i in Placed) {
    ok.chk <- isPFM & (birds$ID == birds$ID[i]);
    Cint   <- c(Cint, diff(age[ok.chk])[-1]);
    ok.src <- isFT  & (birds$ID == birds$ID[i]);
    Sint   <- c(Sint, diff(age[ok.src])[-1]);
  }       # FT  Search Intervals & PFM Check Intervals: cnt, avg, sd
  Ik <- c(n=length(Sint), mu=mean(Sint), sd=sd(Sint));
  Sk <- c(n=length(Cint), mu=mean(Cint), sd=sd(Cint));
  
  #Search dataset
  srch <- data.frame(Id=as.character(birds$ID), Date=Date, Days=nDays,
              Found=birds$Found)[FT,];
  
  if(input=="csv"){
    fn=fname
  }else {
    fn = deparse(substitute(fname))
  }
  invisible(list(scav=scav, srch=srch, Ik=Ik, Sk=Sk,
            NP.Spec=NP.Spec, NP.ID=NP.ID, fn=fn, Info=getId()));
}