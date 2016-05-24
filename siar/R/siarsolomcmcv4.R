siarsolomcmcv4 <-
function(data,sources,corrections=0,concdep=0,iterations=200000,burnin=50000,howmany=10000,thinby=15,prior=rep(1,nrow(sources)),siardata=list(SHOULDRUN=FALSE))
{

    if(siardata$SHOULDRUN==FALSE) {
        siardata <- list()
        siardata$iterations <- iterations
        siardata$burnin <- burnin
        siardata$howmany <- howmany
        siardata$thinby <- thinby
        siardata$TITLE <- "SIAR data"
    }

    # data should be a matrix of two or three columns and numdata rows.
    # If it's 3 columns, the first must be the group number.

    # sources should be a matrix of numiso+1 columns and numsources rows.
    # The first column of sources should be the source names.

    # corrections should be a matrix of 2 columns and numiso rows where the
    # first column is the mean and the second is the sd of the fractionation
    # correction.

    numsources <- nrow(sources)
    numdata <- nrow(data)
    numiso <- (ncol(sources)-1)/2

    if(ncol(data)==numiso+1) {
        #data2 <- data[, 2:(numiso+1)]
	      #data2 <- data[,2:3]
	      data2 <- data.matrix(data[, 2:(numiso + 1)])
        numgroups <- max(data[,1])
        startgroup <- as.vector(c(0,cumsum(table(data[,1])))+1)[1:numgroups]
        endgroup <- as.vector(cumsum(table(data[,1])))
    } else {
        numgroups <- 1
        data2 <- data
        startgroup <- 1
        endgroup <- numdata
    }

    sourcenames <- as.character(sources[,1])
    sourcedata <- sources[,2:(2*numiso+1)]
    parameters <- matrix(1,ncol=(numsources+numiso)*numgroups,nrow=(siardata$iterations-siardata$burnin)/siardata$thinby)

    if(!is.data.frame(corrections)) {
        correctionsdata <- matrix(0,ncol=2*numiso,nrow=numsources)
    } else {
        correctionsdata <- matrix(as.double(as.matrix(corrections[,2:(2*numiso+1)])),nrow=numsources)
    }
    if(!is.data.frame(concdep)) {
        concdepdata <- matrix(1,ncol=2*numiso,nrow=numsources) 
    } else {
        concdepdata <- matrix(as.double(as.matrix(concdep[,2:(2*numiso+1)])),nrow=numsources)
    }
    
    BAD <- FALSE
    if(round((siardata$iterations-siardata$burnin)/siardata$thinby)!=(siardata$iterations-siardata$burnin)/siardata$thinby) {
        cat("Error in iterations, burnin or thinby: (iterations-burnin)/thinby must be an integer. \n \n")
        BAD <- TRUE
    }
    if(!is.numeric(data2[,1])) {
        cat("Error in the target file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if(!is.numeric(sourcedata[,1])) {
        cat("Error in the sources file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if(!is.numeric(correctionsdata[,1])) {
        cat("Error in the corrections file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if(!is.numeric(concdepdata[,1])) {
        cat("Error in the concentration dependence file - check this is numeric. \n \n")
        BAD <- TRUE
    }
    if(numgroups>32) {
      cat("Error: too many groups - take some of the groups out before running. \n")
      BAD <- TRUE
    }

    # Finally, sort out the deterministic problems from the others
    #for(i in 1:numgroups) {
    #    for(j in startgroup[i]:endgroup[i]) {
    #        if(numsources*(endgroup[i]-startgroup[i]+1)==numiso) {
    #            if(numgroups>1) {
    #                cat(paste("Group",i,": this is a deterministic problem and thus not suitable for siar. \n"))
    #            } else {
    #                cat("This is a deterministic problem and thus not suitable for siar. \n")
    #            }
    #            BAD <- TRUE
    #        }
            #if(numsources*(endgroup[i]-startgroup[i]+1)<numiso) {
            #    if(numgroups > 1) {
            #        cat(paste("Group",i,": this is an insoluble problem and thus not suitable for siar. \n"))
            #        cat("The number of isotopes is less than the number of targets times the number of sources. \n")
            #    } else {
            #        cat(paste("This is an insoluble problem and thus not suitable for siar. \n"))
            #        cat("The number of isotopes is less than the number of targets times the number of sources. \n")
            #    }
            #    BAD <- TRUE
            #}
    #    }
    #}

    if(BAD==FALSE) {
        tempout <- .C("siarsolomcmcv4",as.integer(numdata),as.integer(numsources),
                      as.integer(numiso),as.integer(numgroups),as.integer(startgroup),
                      as.integer(endgroup),as.integer(siardata$iterations),
                      as.integer(siardata$burnin),as.integer(siardata$howmany),
                      as.integer(siardata$thinby),as.double(prior),as.double(data2),
                      as.double(concdepdata),as.double(as.matrix(sourcedata)),
                      as.double(correctionsdata),as.double(parameters))
        #,PACKAGE="siar")
        parameters <- matrix(tempout[[16]], ncol = (numsources + numiso) * numgroups,
                             nrow = (siardata$iterations - siardata$burnin)/siardata$thinby)
        if(numgroups==1) {
            colnames(parameters) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
        } else {
            colnames(parameters) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,numiso),sep="")),"G",sep=""),times=numgroups),sort(rep(seq(1,numgroups),times=numsources+numiso)),sep="")
        }

        # Get rid of unwanted columns (set to zero)
        #if(numgroups==1) {
        #  newtempout <- newtempout[,1:numsources]
        #} else {
          # probably a quicker way of doing this but can't think of it right now.
        #  goodcols <- 1:numsources
        #  for(j in 1:(numgroups-1)) {
        #    goodcols <- c(goodcols,(j*(numsources+numiso)+1):(j*(numsources+numiso)+numsources))
        #  }
        #  newtempout <- newtempout[,goodcols]
        #}

        return(list(targets=data,sources=sources,corrections=corrections,concdep=concdep,PATH=siardata$PATH,TITLE=siardata$TITLE,numgroups=tempout[[4]],numdata=tempout[[1]],numsources=tempout[[2]],numiso=tempout[[3]],SHOULDRUN=TRUE,GRAPHSONLY=FALSE,EXIT=FALSE,SIARSOLO=TRUE,output=parameters))

    } else {
        cat("Problems with inputs: siar has not been run. \n")
        return(siardata)
    }


}
