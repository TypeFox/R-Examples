deter.iso <-
function(diffs,histbreaks=50000,mzmin=0,mzmax=0.5,cutcount=180,plotit=TRUE){

   #############################################################################
   # assemble "iso" ############################################################
   his<-hist(diffs,breaks=histbreaks,plot=plotit,xlim=c(mzmin,mzmax),xlab="m/z");
   if(plotit==TRUE){abline(h=cutcount,col="red")};
   masstol2<-c((his[[1]][2]-his[[1]][1])/2);
   his<-his$mids[his$counts>cutcount & his$mids<=mzmax & his$mids>=mzmin];
   #his<-his[his>(0+masstol2)]; # not to search peak = peak !
   that<-data.frame(as.character(rep(1:length(his))),his,rep(0,length(his)),rep(0,length(his)),
   rep(1,length(his)),rep(0,length(his)),rep(1,length(his)));
   names(that)<-c("name","dmass","dabund","how_often","#atoms","/C","z");
   elements<-as.character(rep(1:length(his)));
   charges<-c(1);
   number<-length(his);
   iso<-list(c(),that,charges,number,elements);
   names(iso)<-c("list of isotopes","list of isotope masses","charges","number of isotope m/z","elements");
   #############################################################################
   cat("\n For using deter.iso-output in function pattern.search:\n");
   cat("(1) set deter=TRUE.\n");
   cat("(2) set ppm=FALSE.\n");
   cat(paste("(3) use mztol=",masstol2,".\n\n",sep=""));
   #############################################################################
   return(iso)
}
