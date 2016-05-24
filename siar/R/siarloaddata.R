siarloaddata <-
function(siarversion) {

choices2 <- c("Load data in from files","Load in R objects","Load in previous output")
title <- "The available options are:"
choose2 <- menu(choices2,title = title)

############################################################################

if(choose2==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))

############################################################################
# Load in data from files

if(choose2==1) {

cat("To run siar, you need to have created at least 2 text files. \n")
cat("The first must contain the target isotope measurements in either \n")
cat("two columns with no group number or 3 columns with a group label. \n")
cat("The second file must contain a column of the different source names \n")
cat("followed by the isotope measurements for each in a seperate columns. \n \n")
cat("Optionally, a third file can be created which contains the fractionation \n")
cat("correction means and standard deviations for each isotope. \n \n")
cat("See the demo and the included data files for more information \n")
cat("on the data input format. \n \n")

BADPATH <- TRUE
while(BADPATH == TRUE) {
    cat("First input the directory at which the files can be found: \n")
    PATH <- scan(what="",nlines=1,quiet=TRUE)
    while(length(PATH)==0) PATH <- scan(what="",nlines=1,quiet=TRUE)
    if(PATH==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(PATH)) {
        BADPATH <- FALSE
    } else {
        cat("Cannot find this directory. Check your typing. \n")
    }
}

BADDATA <- TRUE
while(BADDATA == TRUE) {
    cat("Now input the name of the target isotope file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(DATAFILE)==0) DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(DATAFILE==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(paste(PATH,"/",DATAFILE,sep=""))) {
        BADDATA <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }
}

BADSOURCES <- TRUE
while(BADSOURCES == TRUE) {
    cat("Now input the name of the source isotope file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")    
    SOURCEFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(SOURCEFILE)==0) SOURCEFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(SOURCEFILE==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(paste(PATH,"/",SOURCEFILE,sep=""))) {
        BADSOURCES <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }
}

BADCORRECTIONS <- TRUE
while(BADCORRECTIONS == TRUE) {
    cat("Now input the name of the fractionation corrections file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    cat("or leave blank to use pre-corrected values \n")
    CORRECTIONSFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(length(CORRECTIONSFILE)==0) CORRECTIONSFILE <- -999
    if(CORRECTIONSFILE==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(paste(PATH,"/",CORRECTIONSFILE,sep=""))) {
        BADCORRECTIONS <- FALSE
    } else {
        if(CORRECTIONSFILE == -999) {
            BADCORRECTIONS <- FALSE
            corrections <- matrix(0,nrow=1,ncol=1)
        }
        if(CORRECTIONSFILE !=-999) cat("Cannot find this file, check your typing \n")
    }
}

BADCONCDEP <- TRUE
while(BADCONCDEP == TRUE) {
    cat("Now input the name of the concentration dependence file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    cat("or leave blank to ignore concentration dependence. \n")
    CONCDEPFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(length(CONCDEPFILE)==0) CONCDEPFILE <- -999
    if(CONCDEPFILE==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(paste(PATH,"/",CONCDEPFILE,sep=""))) {
        BADCONCDEP <- FALSE
    } else {
        if(CONCDEPFILE == -999) {
            BADCONCDEP <- FALSE
            concdep <- matrix(0,nrow=1,ncol=1)
        }
        if(CONCDEPFILE !=-999) cat("Cannot find this file, check your typing \n")
    }
}


cat("Now loading in data... \n")
targets <- as.data.frame(read.table(paste(PATH,"/",DATAFILE,sep=""),header=TRUE))
sources <- as.data.frame(read.table(paste(PATH,"/",SOURCEFILE,sep=""),header=TRUE))
if(CORRECTIONSFILE != -999) corrections <- as.data.frame(read.table(paste(PATH,"/",CORRECTIONSFILE,sep=""),header=TRUE))
if(CONCDEPFILE != -999) concdep <- as.data.frame(read.table(paste(PATH,"/",CONCDEPFILE,sep=""),header=TRUE))
cat("Done \n \n")

# Finally sort everything out so its in proper siar format
numgroups <- 1
if(targets[1,1]%%1 == 0) numgroups <- max(targets[,1])
numsources <-nrow(sources)
numdata <- nrow(targets)
numiso <- (ncol(sources)-1)/2
if(corrections[1,1] == 0) corrections <- matrix(0,nrow=nrow(sources),ncol=2*numiso+1)
if(concdep[1,1] == 0) concdep <- matrix(0,nrow=nrow(sources),ncol=2*numiso+1)
SHOULDRUN <- TRUE
GRAPHSONLY <- FALSE
EXIT <- FALSE
output <- NULL

}

############################################################################
# Load in R objects

if(choose2==2) {

cat("Please enter the name of the object which contains the target data. \n")
cat("\n")
dataexists <- FALSE
while(dataexists == FALSE) {
    datatemp <- scan(what="",nlines=1,quiet=TRUE)
    while(length(datatemp)==0) datatemp <- scan(what="",nlines=1,quiet=TRUE)
    if(datatemp==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(!exists(datatemp)) {
        cat("Object not found. Try again or Esc to quit. \n")
    } else {
        targets <- get(datatemp)
        dataexists <- TRUE        
    }
}

cat("Now please enter the name of the object which contains the source \n")
cat("isotope details. The first column should be the source names \n")
cat("\n")
sourcesexists <- FALSE
while(sourcesexists == FALSE) {
    sourcestemp <- scan(what="",nlines=1,quiet=TRUE)
    while(length(sourcestemp)==0) sourcestemp <- scan(what="",nlines=1,quiet=TRUE)
    if(sourcestemp==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(!exists(sourcestemp)) {
        cat("Object not found. Try again or Esc to quit. \n")
    } else {
        sources <- get(sourcestemp)
        sourcesexists <- TRUE        
    }
}

cat("Now please enter the name of the object which contains the isotopic \n")
cat("correction mean and standard deviation. Note: if the data are \n")
cat("pre-corrected please leave blank \n")
cat("\n")
correctionsexists <- FALSE
while(correctionsexists == FALSE) {
    correctionstemp <- scan(what="",nlines=1,quiet=TRUE)
    if(length(correctionstemp)==0) {
        corrections <- matrix(0,nrow=1,ncol=1)
        correctionsexists <- TRUE
    } else {
        if(correctionstemp==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
        if(!exists(correctionstemp)) {
            cat("Object not found. Try again or Esc to quit. \n")
        } else {
            corrections <- get(correctionstemp)
            correctionsexists <- TRUE        
        }
    }
}

cat("Now please enter the name of the object which contains the concentration \n")
cat("dependence means and standard deviations. Note: concentration dependence \n")
cat("standard deviations are currently not supported. Leave blank for no \n")
cat("concentration depdendence. \n")
cat("\n")
concdepexists <- FALSE
while(concdepexists == FALSE) {
    concdeptemp <- scan(what="",nlines=1,quiet=TRUE)
    if(length(concdeptemp)==0) {
        concdep <- matrix(0,nrow=1,ncol=1)
        concdepexists <- TRUE
    } else {
        if(concdeptemp==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
        if(!exists(concdeptemp)) {
            cat("Object not found. Try again or Esc to quit. \n")
        } else {
            concdep <- get(concdeptemp)
            concdepexists <- TRUE        
        }
    }
}

# Finally sort everything out so its in proper siar format
numgroups <- 1
if(targets[1,1]%%1 == 0) numgroups <- max(targets[,1])
numsources <-nrow(sources)
numdata <- nrow(targets)
numiso <- (ncol(sources)-1)/2
if(corrections[1,1] == 0) corrections <- matrix(0,nrow=numsources,ncol=2*numiso+1)
PATH <- NULL
SHOULDRUN <- TRUE
GRAPHSONLY <- FALSE
EXIT <- FALSE
output <- NULL

}

############################################################################
# Load in a previous run

if(choose2==3) {

cat("This option allows you to load in the parameters of a previously saved run \n")
cat("so that you can produce graphs, etc. \n \n")
cat("More complex analysis can be run by simply loading this file \n")
cat("into R itself via load(data) (rather than through this menu system).\n \n")

BADOUTPUT <- TRUE
GRAPHSONLY <- TRUE
while(BADOUTPUT == TRUE) {

    cat("Now input the name of the output file including the file extension \n")
    cat("and the directory where it is located (eg","c:\\siar\\data\\output.Rdata)",". \n")

    OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(OUTPUTFILE)==0) OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(OUTPUTFILE==0) return(list(EXIT=FALSE,SHOULDRUN=FALSE))
    if(file.exists(OUTPUTFILE)) {
        BADOUTPUT <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }

}

#output <- as.matrix(read.table(file=OUTPUTFILE,header=TRUE))
load(file=OUTPUTFILE)

# Finally sort everything out so its in proper siar format
siardata <- NULL
targets <- siardata$targets
sources <- siardata$sources
corrections <- siardata$corrections
if(exists(siardata$concdep)) {
    concdep <- siardata$concdep
} else {
    concdep <- matrix(1,nrow=1,ncol=1)
}
PATH <- siardata$PATH
numdata <- siardata$numdata
SHOULDRUN <- siardata$SHOULDRUN
GRAPHSONLY <- siardata$GRAPHSONLY
EXIT <- siardata$EXIT
output <- siardata$output
TITLE <- siardata$TITLE
numgroups <- siardata$numgroups
numdata <- siardata$numdata
numsources <- siardata$numsources
numiso <- siardata$numiso

}

############################################################################
# Enter a title

if(choose2==1 || choose2==2) {
    cat("\n Please enter a name for the data set to be used in plots (or leave blank for default titles). \n")
    TITLE <- scan(what="",nlines=1,quiet=TRUE,sep="\t")
    if(length(TITLE)==0) {
        TITLE <- "SIAR data"
    } else {
        if(TITLE==0) {
          cat("Data not loaded. \n")
          return(EXIT=FALSE,SHOULDRUN=FALSE)
        }
    }
}

return(list(targets=targets,sources=sources,corrections=corrections,concdep=concdep,PATH=PATH,TITLE=TITLE,numgroups=numgroups,numdata=numdata,numsources=numsources,numiso=numiso,SHOULDRUN=SHOULDRUN,GRAPHSONLY=GRAPHSONLY,EXIT=EXIT,output=output))

}
