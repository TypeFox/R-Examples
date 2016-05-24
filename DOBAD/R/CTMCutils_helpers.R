

#graph a CTMC
#CTMC has $states, $times, $T
#save for latex'ing via trellis graphics, if a filename is passed.
#filename should not have an extension nor a period '.'.
#a '.ps' and a '.pdf' will be created.
graph.CTMC <- function(CTMC, filename=NA, height=6, width=4.5,
                       xlab="time",ylab="State",...){
  if (inherits(CTMC, "CTMC")){  #better way to do this?
    times <- getTimes(CTMC);
    states <- getStates(CTMC);
    T <- getT(CTMC);
  }
  else if (class(CTMC)[1] == "list"){
    times <- CTMC$times;
    states <- CTMC$states;
    T <- CTMC$T;
  }
  n <- length(times)
  times <- c(times,T)
  states <- c(states, states[n])
  plot(times, states, type="s",
       ylab=ylab, xlab=xlab,...);
  if ( is.character(filename)){
    ## require(lattice);## alrady required in package
    trellis.device(postscript, file=paste(filename, ".ps", sep=""),
                   horiz=FALSE,
                   height = height, width=width, title=filename);
    plot(times, states, type="s",
         ylab=ylab, xlab=xlab,...);
    dev.off();
    trellis.device(pdf,  file=paste(filename, ".pdf", sep=""),
                   height = height, width=width, title=filename);
    plot(times, states, type="s",
         ylab=ylab, xlab=xlab,...);
    dev.off();
  }
}


#graph partially observed CTMC.
#same as graph.CTMC except different default for type of plot
graph.CTMC.PO <- function(CTMC, filename=NA,height=6, width=4.5,
                          type="l",...){
#  par(pointsize=.1) #remember to reset at end
  if (inherits(CTMC, "CTMC_PO_1")){
    times <- getTimes(CTMC);
    states <- getStates(CTMC);
    T <- getT(CTMC);
  }
  else if (class(CTMC)[1] == "list"){
    times <- CTMC$times;
    states <- CTMC$states;
    T <- CTMC$T;
  }  
  plot(times, states, type=type,
       ylab="State", xlab="time",
       ...);
  if ( is.character(filename)){
#    filename <- paste(filename, "PO", sep="")
    ##require(lattice);
    trellis.device(postscript, file=filename, horiz=FALSE,
                   height = height, width=width, title=filename);
      plot(times, states, type=type,
       ylab="State", xlab="time",
       ...);
    dev.off();
  }
}

###Use "unlist" rathre than simplify !! ! !   

#takes a list, each of whose entry is a numeric vec length 1 and turns it into a vector
simplify <- function(simpleList){
  tester <- function(x){
    if(length(x)>1 || (mode(x) !="numeric")) {
      print("error");
      stop();
    }
  x[1];
  }
  sapply(simpleList,tester, simplify=TRUE);
}
