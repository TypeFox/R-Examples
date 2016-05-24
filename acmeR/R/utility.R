#'@import foreign
#'@import utils
# -----------------------Utility Function---------------------- #
#       CalWEA2acme() Make ACME dataset from CAL/Wea spreadsheet
#       class2spec()  Find species list for "class" short-cuts like "BATS"
#       spec2class()  Find class for specified species
#       spec2name()   Find common or scientific name of species (code)
#       subspec()     Subset of data for specified species (one or many)
#       ages()        Find carcass ages (in days) at times of search
#       str2str()     Format string to specified width, left justified
#       int2str()     Format integer as string of specified width, right j
#       bleed()       Numbers of observed bleed-through events
#       getId()       Reports info about invocation (date, machine, etc)

##################################################################
svnid <- "$Id: acme.R 128 2015-05-20 10:47:48Z rlw $";
# CalWEA2acme(): Read original CAL/Wea dataset and construct acme
#              spreadsheet from it for use by acme package.
#
CalWEA2acme <- function(fname="acme.csv", tz="PST8PDT") {
  birds   <- read.csv("FinalFullDataSet.csv");
  nrecs   <- dim(birds)[1];
#
  Date    <- as.character(birds$date);
  Time    <- rep("12:00:00 PM",nrecs);  # Not included in Cal-WEA file
#
  ID      <- as.character(birds$AssignedID);
#
  Species <- as.character(birds$Species);
#
  PFM     <- birds$Person=="PFM";
  Placed  <- (birds$DetectionStatus=="P");
#
  Event                <- rep(NA,nrecs);
  Event[Placed]        <- "Place";
  Event[PFM & !Placed] <- "Check";
  Event[!PFM]          <- "Search";
#
  Found                <- (birds$searcherfind > 0);
  Found[is.na(Found)]  <- FALSE;
#
  rv <- data.frame(Date=Date,Time=Time,ID=ID,Species=Species,
                   Event=Event,Found=Found);
  write.table(rv, file=fname, sep=",", qmethod="double", row.names=F);
  invisible(rv);
}

##################################################################
  acme.bats  <- c("BBBA", "LBBA", "SHBA", "UNBA", "UNMY", "UNPI");
  acme.small <- c("ACWO", "AMGO", "AMKE", "AMRO", "BHCO", "BHGR",
                  "BOWA", "BRBL", "BUSH", "EUST", "FOSP", "GCSP",
                  "HOLA", "LEGO", "MODO", "NOFL", "WCSP", "WEME",
                  "WESJ");
  acme.large <- c("AMCR", "BNOW", "CAGU", "CORA", "GHOW", "GOEA",
                  "GREG", "MALL", "ROPI", "RTHA", "SWHA", "TUVU",
                  "UNRA");
  acme.huge  <- c("CAGO", "WITU", "UNDG");
##################################################################
class2spec <- function(spec) {
  # class2spec(): Find species list for "class" short-cuts like "BATS"
  
  if(tolower(spec[1]) == "bats")                   return(acme.bats)
  else if(tolower(spec[1]) == "small")             return(acme.small)
  else if(tolower(spec[1]) %in% c("big","large"))  return(acme.large)
  else if(tolower(spec[1]) %in% c("extra", "huge"))return(acme.huge)
  else if(tolower(spec[1]) == "birds")
         return(c(acme.small, acme.large, acme.huge))
  else if(tolower(spec[1]) == "all")
         return(c(acme.bats, acme.small, acme.large, acme.huge))
  else   return(spec);
}

##################################################################
#
spec2class <- function(cls) { 
  # Inverse of "class2spec": report class for specific species
  
  if(toupper(cls[1]) %in% c("BATS", acme.bats))    return("Bats")
  else if(toupper(cls[1]) %in% c("BIG","LARGE", acme.large))
    return("Big Birds")
  else if(toupper(cls[1]) %in% c("SMALL", acme.small))
    return("Small Birds")
  else if(toupper(cls[1]) %in% c("HUGE", "EXTRA", acme.huge))
    return("Extra-large Birds")
  else
    return(cls[1]);
}

##################################################################

#
spec2name <- function(spec, sci=FALSE, fname="species.csv")
  # Return common or scientific name of species
  
 {
  if(!length(spec)) return("All");
  spec.names <- species.file; #species.file is an internal package dataset
  dbfs <- list.files(pattern="LIST.*DBF");
  if(ld <- length(dbfs)) {
    aou <- read.dbf(dbfs[ld], as.is=TRUE);
    aou.names <- cbind(Class="BIRD", Code=aou$SPEC,
       Sci=aou$SCINAME, Com=aou$COMMONNAME);
    spec.names <- rbind(spec.names, aou.names);
  }
  found      <- spec.names$Code %in% spec;
  if(!any(found)) return(paste(spec,collapse=", "));
  if(sum(found)>1) {
    found.codes <- cbind(which(found),spec.names$Code[found]);
    for(i in 2:sum(found)) {  # Remove duplicates:
      if(found.codes[i,2] %in% found.codes[1:(i-1),2]) {
        found[as.numeric(found.codes[i,1])] <- FALSE;
      }
    }
  }
  if(sci) rv <- unique(spec.names$Sci[found])
  else    rv <- unique(spec.names$Com[found]);
  return(paste(rv, collapse=", "));
}

##################################################################
#
subspec <- function(rd, spec) {
  # Subset for specific species
  
  if(any(tolower(spec) %in% c("all", ""))) return(rd);
  ok     <- which(rd$scav$Species %in% class2spec(spec));
  if(!any(ok)) stop("Sorry, no species \"", spec, "\" found.");
  scav   <- rd$scav[ok,];
  srch   <- rd$srch[rd$srch$Id %in% scav$Id,];
  rd$scav <- scav;
  rd$srch <- srch;           # No Sk
  return(list(scav=scav, srch=srch, Ik=srch.interval(rd),
         Info=rd$Info));
}

##################################################################

#
srch.interval <- function(rd) {
  # Sample mean & variance for lengths of FT Search intervals Ik
  # We skip first search for each carcass, to look at *intervals*.
  
  scav  <- rd$scav;
  srch  <- rd$srch;
  tot   <- dim(scav)[1];
  Ik    <- c(n=0, mu=0, s2=0);              # FT  Search Intervals
  for(i in 1:tot) {
    ok  <- srch$Id==scav$Id[i];             # This carcass?? (T/F)
    m.i <- sum(ok)-1;                       # This carcass search count
    if(m.i > 0) {                           # At least TWO searches...
      age.i    <- diff(as.numeric(difftime(srch$Date[ok], # Carcass ages
                       scav$Placed[i], units="days")));   # at search (days)
      Ik["n"]  <- Ik["n"]  + m.i;           # Increment total search count
      f.i      <- m.i/Ik["n"];              # Fraction for this carcass
      mu.i     <- mean(age.i);              # Average length this carc
      s2.i     <- sum((age.i-mu.i)^2)/m.i;  # MLE for variance, this sample
      delta    <- mu.i - Ik["mu"];          # Innovation
      Ik["mu"] <- Ik["mu"] + f.i * delta;   # Recursive update of mean
      Ik["s2"] <- (1-f.i) * Ik["s2"] +      # Recursive update of variance
                  f.i * s2.i +
                  f.i * (1-f.i) * delta^2;
    }
  }
  return(Ik);
}

##################################################################
#
ages <- function(rd) {
  # ages():  Find ages (in days) at times of search
  
  scav  <- rd$scav;
  srch  <- rd$srch;
  days  <- numeric(dim(srch)[1]);
  ncarc <- dim(scav)[1];
  for(i in 1:ncarc) {
    ok       <- srch$Id==scav$Id[i];         # Entries for this carcass
    days[ok] <- as.numeric(difftime(srch$Date[ok],
                scav$Placed[i], units="days"));
  }
  return(days);
}


##################################################################
# str2str():  Format string x to specified width w
#
str2str <- function(x, w) {
  substring(paste(x,paste(rep(' ',w),collapse=''),collapse=''),1,w);
}
##################################################################
# int2str():  Format integer x to specified width w
#
int2str <- function(x, w) {
  rv <- paste(paste(rep(' ',w),collapse=''),x,collapse='');
  len <- nchar(rv); return(substring(rv,len+1-w,len));
}

##################################################################

#
bleed <- function(rd) {
  # Find observed bleed-through counts for each carcass,
  # i.e., number of discoveries after a search failure
  
  scav   <- rd$scav;
  srch   <- rd$srch;
  n.carc <- dim(scav)[1];
  rv     <- numeric(n.carc);        # Number of succ after a miss
  for(i in 1:n.carc) {
    ok    <- srch$Id==scav$Id[i];   # Entries for this carcass
    fnd   <- srch$Found[ok];
    rv[i] <- sum(fnd[-(1:min(sum(ok),which(!fnd)))]); 
  }
  return(rv);
}

#######################################################################
#
getId <- function() {
  # Gather info about user, system, date, environment, etc. to save with
  # program output, making it easier to keep track of ouptut files.
  # Try to isolate Windows- and Unix- dependencies.
  
  acme  <- paste("ACME build ", strsplit(svnid," ")[[1]][3], ": ",
                strsplit(svnid," ")[[1]][4], sep="");
  when  <- Sys.time();
  info  <- Sys.info();
  what  <- info[1];   # sysname
  where <- info[4];   # nodename
  who   <- info[7];
  Rver  <- version[13][[1]];
  call  <- sys.call(sys.parent(sys.parent()));
  switch(what,
         "Windows" = { os <- system("cmd /c ver",intern=TRUE,invisible=TRUE)[2]; },
         "Linux"   = { os <- system("uname -rps",intern=TRUE); },
         "Darwin"  = { os <- paste("OSX:",R.version$platform); },
         os <- "Unknown");
  
  return(list("acme"=acme, "Rver"=Rver, "os"=os, "mach"=where,
              "arch"=what, "who"=who, "when"=when, "call"=call));
}