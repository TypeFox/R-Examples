siarsaveoutput <-
function(siardata) {

if(siardata$SHOULDRUN==FALSE) {
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

BADFILE <- TRUE
while(BADFILE == TRUE) {
    cat("This option will save all the model details to a .Rdata file.\n")
    cat("It can be loaded back in via siarmenu() or by the command loadsiardata(file).\n")
    cat("Enter a directory location where the output parameters will reside: \n")

    outputfileloc <- scan(what="",nlines=1,quiet=TRUE)
    while(length(outputfileloc)==0) outputfileloc <- scan(what="",nlines=1,quiet=TRUE)
    if(outputfileloc==0) return(NULL)

    if(!file.exists(outputfileloc)) {
        cat("This location doesn't exist, check your typing \n")
    } else {
        cat("Please enter a filename: \n")
        outputfilename <- scan(what="",nlines=1,quiet=TRUE)
        while(length(outputfilename)==0) outputfilename <- scan(what="",nlines=1,quiet=TRUE)
        if(outputfilename==0) return(NULL)

        if(length(siardata$sources)>0) {
            sourcenames <- as.character(siardata$sources[,1])
        } else {
            sourcenames <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,siardata$numiso),sep="")),"G",sep=""),times=siardata$numgroups),sort(rep(seq(1,siardata$numgroups),times=siardata$numsources+siardata$numiso)),sep="")
        }

        # Get some column names
        if(siardata$targets[1,1]%%1!=0) {
            siarcolnames <- c(sourcenames,paste("SD",seq(1,siardata$numiso),sep=""))
        } else {
            siarcolnames <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,siardata$numiso),sep="")),"G",sep=""),times=siardata$numgroups),sort(rep(seq(1,siardata$numgroups),times=siardata$numsources+siardata$numiso)),sep="")     
        }

        
        cat("Writing output ... \n")
        save(siardata,file=paste(outputfileloc,"/",outputfilename,".Rdata",sep=""))
        BADFILE <- FALSE
        cat("Output created. \n \n ")
        
        cat("Press <Enter> to continue")
        readline()
        invisible()

    }

}





}
