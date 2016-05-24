plotadduct <-
function(adduct){

    def.par <- par(no.readonly = TRUE) # save default, for resetting...
    sc<-close.screen();if(sc[1]!=FALSE){for(m in 1:length(sc)){close.screen(sc[m])}};
    split.screen(c(2,1));
    screen(1) ##################################################################
    par(mar=c(0,4,2,2));
    plot.window(xlim=c(0,(length(adduct[[4]][,2])+1)),ylim=c(0,max((adduct[[4]][,2]))));
    axis(2);box();
    title(ylab="Absolute frequency",xlab="Adduct");
    for(i in 1:length(adduct[[4]][,2])){
        points(i,adduct[[4]][i,2],type="h",lwd=5);
        points(i,adduct[[4]][i,2],type="h",lwd=5,col="red");
    };
    rm(i);
    screen(2) ##################################################################
    here<-matrix(nrow=length(adduct[[4]][,2]),ncol=max(adduct[[4]][,2]),0);
    count<-rep(1,length(adduct[[4]][,2]))
    for(i in 1:length(adduct[[3]][,1])){
           this1<-as.numeric(strsplit(as.character(adduct[[3]][i,2]),",")[[1]]);
           this2<-strsplit(as.character(adduct[[3]][i,3]),"/")[[1]];
           for(j in 1:length(this1)){
                      here[adduct[[4]][,1]==this2[j],count[adduct[[4]][,1]==this2[j]]]<-log10(adduct[[1]][this1[j],2]);
                      if(is.na(here[adduct[[4]][,1]==this2[j],count[adduct[[4]][,1]==this2[j]]])){stop("NA reached!!!")}
                      count[adduct[[4]][,1]==this2[j]]<-count[adduct[[4]][,1]==this2[j]]+1
           };
    };
    rm(i);
    par(mar=c(4,4,0,2));
    plot.window(xlim=c(0,(length(adduct[[4]][,2])+1)),ylim=c(0,max(here[,])));   # max(here[,])
    axis(2);
    title(ylab="log( peak intensities )",xlab="Adducts");
    for(i in 1:length(adduct[[4]][,2])){
      boxplot((here[i,1:(count[i]-1)]),at=i,add=TRUE)
      mtext(adduct[[4]][i,1],at=i,side=1,cex=0.7,line=1)
    };
    rm(i);
    close.screen(all.screens = TRUE);
    # output ###################################################################
    here<-matrix(nrow=length(adduct[[4]][,2]),ncol=length(adduct[[4]][,2]),0);
    rownames(here)=adduct[[4]][,1];
    colnames(here)=adduct[[4]][,1];
    for(i in 1:length(adduct[[1]][,7])){if(adduct[[1]][i,7]!="none"){
              this1<-strsplit(as.character(adduct[[1]][i,7]),"//")[[1]];
              for(j in 1:length(this1)){
                    this2<-strsplit(as.character(this1[j]),"<->")[[1]];
                    here[adduct[[4]][,1]==this2[1],adduct[[4]][,1]==this2[2]]<-c(here[adduct[[4]][,1]==this2[1],adduct[[4]][,1]==this2[2]]+1);
    }}};
    diag(here)<-adduct[[4]][,2]
    rm(i);
    cat("Hit matrix:\n");
    return(here)
    ############################################################################
    par<-par(def.par)#- reset to default

}
