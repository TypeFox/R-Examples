plotdiff <-
function(
	peaklist,
	histbreaks=10000,
	rttol=c(0,0),
	mztol=c(0,100),
	plotit=TRUE
){

  ##############################################################################
  # screen over all differences ################################################
  cat("\n Screen m/z ... please wait!\n");
  if(length(rttol)!=2){stop("rttol must have a lower and an upper bound!")}
  if(rttol[1]>rttol[2]){stop("minimum > maximum for rttol!")}
  if(mztol[1]>mztol[2]){stop("minimum > maximum for mztol!")}
  if(length(mztol)!=2){stop("masstol must have a lower and an upper bound!")}
  if(!is.data.frame(peaklist)){stop("peaklist must be a data.frame")}
  if(length(peaklist[1,])>3){stop("peaklist with > 3 columns not allowed")}
  if(!is.numeric(peaklist[,1]) || !is.numeric(peaklist[,2]) || !is.numeric(peaklist[,3]) ){stop("peaklist columns not numeric")}
  difs<-c();
  fin<-c(length(peaklist[,1]));
  getit<-seq(1,fin,1);
  for(i in 1:fin){
      getit1<-getit[
              peaklist[,3]>=(peaklist[i,3]+rttol[1])  &
              peaklist[,3]<=(peaklist[i,3]+rttol[2])  &
              peaklist[,1]>=(peaklist[i,1]+mztol[1])  &
              peaklist[,1]<=(peaklist[i,1]+mztol[2]) 
              ]; 
      difs<-c(difs,peaklist[getit1,1]-peaklist[i,1]);
  }
  ##############################################################################
  # plot! ######################################################################
  if(plotit==TRUE){
    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    #plot.new();
    sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
    split.screen(c(3,1));
    screen(1);
    par(mar=c(1,4,1,1));
    plot.new();
    hist(difs,breaks=histbreaks,xlim=c(mztol[1],mztol[2]),xlab="",main="Delta mass distributions");
    screen(2);
    par(mar=c(1,4,1,1));
    plot.new();
    hist(difs,breaks=histbreaks,xlim=c(mztol[1],mztol[2]/2),xlab="",main="");
    screen(3);
    par(mar=c(4,4,1,1));
    plot.new();
    hist(difs,breaks=histbreaks,xlim=c(mztol[1],mztol[2]/20),xlab="Delta m/z",main="");
    close.screen(all.screens = TRUE);
    par<-def.par # - reset to default
  }
  ##############################################################################
  return(difs);

}
