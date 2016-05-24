siarhistograms <-
function(siardata,siarversion=0,legloc='topright') {

if(siardata$SHOULDRUN==FALSE && siardata$GRAPHSONLY==FALSE) {
    cat("You must load in some data first (via option 1) in order to use \n")
    cat("this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(NULL)
}

if(length(siardata$output)==0) {
    cat("No output found - check that you have run the SIAR model. \n \n")
    return(NULL)
}

cat("Plots of single groups proportions. \n")

if(siardata$numgroups>1) {
    cat("Enter the group number of the proportions you wish to plot \n")
    #cat("Click on the graph to position the legend... \n")
    BADGROUP <- TRUE
    while(BADGROUP==TRUE) {
      groupnum <- as.integer(scan(what="",nlines=1,quiet=TRUE))
      if(length(groupnum)>0) {
        BADGROUP <- FALSE
        if(groupnum>siardata$numgroups) {
          BADGROUP <- TRUE
          cat("Group number out of range. \n")
        }
      }
    }
} else {
    groupnum <- 1
}

title <- "Do you require each plot on a seperate graph or all on the same one?"
choices <- c("Each on a seperate graph","All together on one graph")
choose <- menu(choices,title = title)

cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(length(siardata$sources)>0) {
    sourcenames <- as.character(siardata$sources[,1])
} else {
    sourcenames <- strsplit(colnames(siardata$output[,((groupnum-1)*(siardata$numsources+siardata$numiso)+1):(groupnum*(siardata$numsources+siardata$numiso)-siardata$numiso)]),paste("G",groupnum,sep=""))
}

# Get the right dimensions of pars
usepars <- siardata$output[,((groupnum-1)*(siardata$numsources+siardata$numiso)+1):(groupnum*(siardata$numsources+siardata$numiso))]

mybreaks <- seq(0,1,length=50)
halfwidth <- diff(mybreaks)[1]/2
top <- 0
for(j in 1:siardata$numsources) {
    top <- max(c(top,max(hist(usepars[,j],plot=FALSE,breaks=mybreaks)$density)))
}

if(choose==2) {

if(siardata$TITLE!="SIAR data") {
    if(siardata$numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(siardata$TITLE,": proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(siardata$numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(siardata$TITLE,": proportion densities",sep=""),xlab="proportion",ylab="density")
} else {
    if(siardata$numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(siardata$numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main="Proportion densities",xlab="proportion",ylab="density")
}
if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:siardata$numsources) {
    Ans <- hist(usepars[,j],plot=FALSE,breaks=mybreaks)
    for(k in 1:length(Ans$mids)) {
        lines(c(Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(siardata$numsources+1)/2,lend=1)
        lines(c(Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(siardata$numsources+1)/2,lend=1)
    }
}

legend(legloc,legend=sourcenames,col=seq(1,siardata$numsources),lty=1,lwd=3,bty="n")

}

if(choose==1) {
devAskNewPage(ask=TRUE)
  
for(j in 1:siardata$numsources) {
if(siardata$TITLE!="SIAR data") {
    if(siardata$numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(siardata$TITLE,": proportion densities for group ",groupnum,": ",sourcenames[j],sep=""),xlab="proportion",ylab="density")
    if(siardata$numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(siardata$TITLE,": proportion densities: ",sourcenames[j],sep=""),xlab="proportion",ylab="density")
} else {
    if(siardata$numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities for group ",groupnum,": ",sourcenames[j],sep=""),xlab="proportion",ylab="density")
    if(siardata$numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities: ",sourcenames[j],sep=""),xlab="proportion",ylab="density")
}
    if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

    Ans <- hist(usepars[,j],plot=FALSE,breaks=mybreaks)
    for(k in 1:length(Ans$mids)) {
        lines(c(Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(siardata$numsources+1)/2,lend=1)
        lines(c(Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((siardata$numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(siardata$numsources+1)/2,lend=1)
    }
}
}



}
