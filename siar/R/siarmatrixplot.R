siarmatrixplot <-
function(siardata,siarversion=0) {

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

cat("Matrix plot of groups proportions. \n")

if(siardata$numgroups>1) {
    cat("Enter the group number of the proportions you wish to plot \n")
    cat("Click on the graph to position the legend... \n")
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

cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(length(siardata$sources)>0) {
    sourcenames <- as.character(siardata$sources[,1])
} else {
    sourcenames <- strsplit(colnames(siardata$output[,((groupnum-1)*(siardata$numsources+siardata$numiso)+1):(groupnum*(siardata$numsources+siardata$numiso)-siardata$numiso)]),paste("G",groupnum,sep=""))
}


# Get some column names
if(length(siardata$targets)>0) {
    if(siardata$numgroups==1) {
        colnames(siardata$output) <- c(sourcenames,paste("SD",seq(1,siardata$numiso),sep=""))
    } else {
        colnames(siardata$output) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,siardata$numiso),sep="")),"G",sep=""),times=siardata$numgroups),sort(rep(seq(1,siardata$numgroups),times=siardata$numsources+siardata$numiso)),sep="")     
    }
}

# Get the right dimensions of pars
usepars <- siardata$output[,((groupnum-1)*(siardata$numsources+siardata$numiso)+1):(groupnum*(siardata$numsources+siardata$numiso))]

if(length(siardata$TITLE) > 0) {
    if(siardata$TITLE!="SIAR data") {
        if(siardata$numgroups > 1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main=paste(siardata$TITLE,": matrix plot of proportions for group ",groupnum,sep=""),diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
        if(siardata$numgroups ==1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main=paste(siardata$TITLE,": matrix plot of proportions",sep=""),diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
    } else {
        if(siardata$numgroups > 1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main=paste("Matrix plot of proportions for group ",groupnum,sep=""),diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
        if(siardata$numgroups ==1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main="Matrix plot of proportions",diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
    }
} else {
    if(siardata$numgroups > 1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main=paste("Matrix plot of proportions for group ",groupnum,sep=""),diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
    if(siardata$numgroups ==1) pairs(usepars[,1:siardata$numsources],xlim=c(0,1),ylim=c(0,1),main="Matrix plot of proportions",diag.panel=panelhist,lower.panel=panelcor,upper.panel=panelcontour)
}

if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

}
