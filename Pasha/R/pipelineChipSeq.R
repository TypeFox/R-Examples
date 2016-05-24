#############################################
##                                         ##
## Pipeline for chipSeq/RNASeq experiments ##
##                                         ##
#############################################
#
# Author  : Romain Fenouil
# Date    : 2011-06-10
#
# THIS VERSION MUST BE USED WITH R >= 2.13.0  AND SHORTREAD LIBRARY >= 1.10.1
# THERE IS AN INCONSISTENCY PROBLEM WITH THE PRECEDENT VERSIONS OF THE LIBRARY (INVERTED STRANDS WHEN USING READALIGNED FUNCTION)


###################################
## Multiprocessors library choice (deprecated since 'parallel' has been included in R)
#
# After several tries, I chose 'multicore' library. Snow and snowfall are designed for being flexible and useable on networks also. 
# It implies that they don't share memory but need copy of objects to each new thread. Since we have to work on big objects, it's
# better to find a library which allows memory share or at least that limit the number of copies in memory. This is the case of
# 'multicore' which uses 'fork' and is designed to be used on a single computer with multi processors or cores. Since the OS is
# supposed to provide Copy-On-Write mechanism , it should avoid copying these big objects to children processes if we don't modify
# them in the function called. One should care about not modifying the big objects in the parallel functions but also limiting the
# size of returned abjects in order to limit the amount of data to copy throught the pipes between the parent process and the children
# (Here it could be very long). I could have used the library 'foreach' which allows more flexibility since it can use either 'snow'
# or 'multicore' as a backend. But 1 -> 'foreach' doesnt provide a ready to use 'lapply' parallel equivalent function that I needed.
# 2-> One should NEVER use this pipeline with 'snow' or network computing since copy of this big objects would be terribly unefficient !
# That's why multicore is the best suited (not ideal though), I just use mclapply directly in a very naive way, trying to minimize the
# number of objects copied. It is also easy to handle text messages from children using 'silent' parameter to 'mclapply'.

#library(multicore)
#library(ShortRead)


# A function that encapsulates the elements of a list (usually useful for big objects) in environments to give a list of environments (pass by environments through mclapply mechanism)
.encapsulateListElementsInEnv=function(listObject)
{
    #cat("\n Encapsulating list in environments")
    return(lapply(listObject,function(x)
                    {
                        capsule <- new.env() # Create the new environment
                        assign("value", x, envir=capsule) # Put the object in the environment
                    }))
}


# A function that checks for the preexistence of an eventual folder or a file with the desired name before actually creating the folder
.safeCreateFolder=function(folderToCreate)
{
    if(!file.exists(folderToCreate))
    {
        dir.create(folderToCreate, recursive=TRUE)
    }
    else
    {
        if(!file.info(folderToCreate)$isdir)
        {
            stop("Cannot create folder for storing results, a file (not a folder) with the same name already exists :", folderToCreate)
        }
    }
}


# Equivalent of deprecated pileup function from ShortRead using the 'coverage' function from Irange
# Inspired from : https://stat.ethz.ch/pipermail/bioconductor/2009-March/026760.html
.newPileup=function(start, fragLength, dir="+", readLength=fragLength)
{
    x_start <- start
    ii <- dir == "-"  # which ranges need to be shifted
    ii_readLength <- if (length(readLength) > 1) readLength[ii] else readLength
    ii_fragLength <- if (length(fragLength) > 1) fragLength[ii] else fragLength
    x_start[ii] <- x_start[ii] + ii_readLength - ii_fragLength #+ 1L
    
    x <- IRanges(start=x_start, width=fragLength)
    
    return(coverage(x))
}


# It is in charge of removing the NULL list elements or objects of length(0)
# Useful to get error messages in forks while getting results from mclapply
.checkErrorsInFork=function(mclapplyResult)
{

    # Find the try-error type class in the result
    hasFailed <- sapply(mclapplyResult, is, "try-error")
    
    if(any(hasFailed))
    {
        stop(paste("\nAt least one job failed with the folling error message(s) :\n  ",
                        paste(sapply(mclapplyResult[hasFailed], as.character), collapse="\n  "),
                        sep=""))
    }
    
    # The other error finding  mechanism (grep in text message) is kept for compatibility
    isNULL_mclapplyResult <- sapply(mclapplyResult, is.null)
    isEmpty_mclapplyResult <- (sapply(mclapplyResult, length)==0)
    
    mclapplyResult <- mclapplyResult[!(isNULL_mclapplyResult | isEmpty_mclapplyResult)]
    
    if(length(mclapplyResult)>0)
    {
        
        # Check for error messages
        isAtomicCharacterResult <- sapply(mclapplyResult, function(x) {return(is.character(x) && length(x==1))})
        if(any(isAtomicCharacterResult))
        {
            mclapplyResult_Character <- mclapplyResult[isAtomicCharacterResult]
            isErrorMessage <- sapply(mclapplyResult_Character, function(x) {return(grepl("error",x))}) 
            
            if(any(isErrorMessage))
            {
                
                cat("\nError(s) in forks :\n")
                sapply(mclapplyResult_Character[isErrorMessage], cat, "\n")
                stop("\nErrors were detected in at least one thread...")
            }
        }
        return(mclapplyResult)
    }
    else
    {
        stop("The resulting list from last loop is empty")
    }
}


# Small function used to display the execution time of the analysis steps
.printTimes=function(returnList)
{
    mapply(function(x, expName)
            {
                cat("\n\n---------------") 
                cat("\nEXP :",expName) 
                mapply(function(x, operationName)
                        {
                            cat("\n", operationName, ":", x["elapsed"], "sec(s) -->", x["elapsed"]/60, "min(s)")
                        }
                        ,x[["execTime"]], names(x[["execTime"]]))
            },
            returnList, names(returnList))
    cat("\n")
    return(NULL)
}


# This is a more generic function that will allow the user to define its own validity function
.checkComplexParam=function(param, paramName, expList, validValuesFct=function(...){return(TRUE)})
{
    if(!is.function(validValuesFct)) stop("validValuesFct must be a a function that tells whether the parameter values are valid or not")
    
    # these param must be atomic, vector, or  list (with same names as experiments names)
    if(is.list(param))
    {
        if(!all(names(expList) %in% names(param))) 
        {
            stop(paste(paramName," parameter is a list but the elements names don't match the experiments names",sep=""))
        }
        else 
        {
            if(!all(sapply(param, validValuesFct), na.rm=TRUE))
            {
                stop(paste(paramName," is a correct list but some elements contain invalid values",sep=""))
            }
            
            # param was already well formatted and the parameters fit the expectations specified by validValuesFct, return it
            return(param)
            
        }
    }
    else
    {
        if(!validValuesFct(param)) 
        {
            stop(paste(paramName,"parameter values are not valid"))
        }
        
        # If it's a vector, we create the list repeating these parameters list for all experiments
        return(lapply(expList, function(x, values){return(values)}, param))
        
    }
}


# addalpha() - adapted from "https://github.com/mylesmharrison/colorRampPaletteAlpha/blob/master/colorRampPaletteAlpha.R"
.addAlpha <- function(colors, alpha=1.0) 
{
    r <- col2rgb(colors, alpha=TRUE)
    # Apply alpha
    r <- r/255
    r[4,] <- alpha
    return(rgb(r[1,], r[2,], r[3,], r[4,]))
}


# Function that plot a density and highligth specific ranges in the distribution (separate by a vertical line and/or fill the range under the curve)
.plotDensityRanges=function(x, ranges=NULL, ranges.col="gray", ranges.alpha=0.3, labels.col="black", include.lower=FALSE, include.upper=TRUE, sep.col="black", sep.lty=2, sep.lwd=1, ...)
{
    
    # Compute and plot the actual density function
    densityPoints <- density(x, n=2048)
    plot(densityPoints, ...)
    
    if(!is.null(ranges))
    {
        
        # Convert 'ranges' to a list to allow the same processing for single and multi ranges
        if(!class(ranges)=="list") ranges=list(ranges)
        
        
        # Summarize all parameters in a data frame to take advantage of recycling mechanism
        parametersDF <- as.data.frame(cbind(ranges = ranges,
                        ranges.col    = .addAlpha(ranges.col, ranges.alpha),
                        labels.col    = labels.col,
                        include.lower = include.lower,
                        include.upper = include.upper))
        
        
        # Keep track of the vertical separation drawn in order to avoid plotting twice on the same coordinate.
        alreadyDrawnSep <- NULL
        
        # Plot all ranges
        for(index in 1:nrow(parametersDF))
        {
            # Get the parameters for the current range in the dataframe
            currentRange      <-  parametersDF[["ranges"]][[index]]
            ranges.col        <-  parametersDF[["ranges.col"]][[index]]
            labels.col        <-  parametersDF[["labels.col"]][[index]]
            include.lower     <-  parametersDF[["include.lower"]][[index]]
            include.upper     <-  parametersDF[["include.upper"]][[index]]
            
            
            # Get the curves coordinates with the closest values from the desired range
            leftInd <- which.min(abs(min(currentRange)-densityPoints$x))
            rightInd <- which.min(abs(max(currentRange)-densityPoints$x))
            
            if(leftInd!=rightInd)
            {
                # Get the coordinates of the range limits (left and right vertical lines)
                coordsRangeLim.x <- c(densityPoints$x[leftInd],densityPoints$x[rightInd])
                coordsRangeLim.y <- c(densityPoints$y[leftInd],densityPoints$y[rightInd])
                
                # Avoid to plot several time the separation lines (detect the ones already plotted)
                alreadyPlotted <- (coordsRangeLim.x %in% alreadyDrawnSep)
                alreadyDrawnSep <- c(alreadyDrawnSep,coordsRangeLim.x[!alreadyPlotted])
                
                # Plot vertical line of separation
                abline(v=coordsRangeLim.x[!alreadyPlotted], col=sep.col, lty=sep.lty, lwd=sep.lwd)
                
                # Define the range name
                labelText <- paste(if(include.lower) "[" else '(',min(currentRange), ",", max(currentRange), if(include.upper) "]" else ')', sep="")
                
                # Plot the label at the middle of the range and defined y coordinate and relative position
                text(x=densityPoints$x[mean(c(leftInd,rightInd))], y=0, labels=paste("    ", labelText, sep=""), adj=0, srt=90, col=labels.col, cex=0.5)
                
                
                # Get all the points corresponding to this part of the curve
                densityPointsInRange.x <- densityPoints$x[leftInd:rightInd]
                densityPointsInRange.y <- densityPoints$y[leftInd:rightInd]
                
                # Add begin and end points at bottom
                densityPointsInRange.x <- c(densityPointsInRange.x[1], densityPointsInRange.x, densityPointsInRange.x[length(densityPointsInRange.x)])
                densityPointsInRange.y <- c(0,densityPointsInRange.y,0)
                polygon(x=densityPointsInRange.x, y=densityPointsInRange.y, col=ranges.col, border=NA)
            }
        }
    }
    
    return(densityPoints)
}


# This function plots a double symmetric histogram
.plotBackToBackBarplot=function(top, bottom, topLab="", bottomLab="" ,main="", cutLineColumn=NULL, gridCol=NULL, outsideLayout=FALSE, plotLegend=TRUE,...)
{
    if(!outsideLayout)
    {
        layout(1:2, heights=c(1,1))
    }
    
    if(is.null(top) || (length(top)==0))
    {
        frame()
    }
    else
    {
        par(mar=c(1.5,4.5,3,3))
        barplot(top, main=main, axes=FALSE, legend.text=plotLegend, ylab=topLab,...)
        if(is.numeric(cutLineColumn)) abline(v=cutLineColumn*1.2+0.1, lwd=2, lty=2, col="red")
        axis(2)
        grid(col=gridCol, nx=0, ny=NULL)
    }
    
    
    
    if(is.null(bottom) || (length(bottom)==0))
    {
        frame()
    }
    else
    {
        par(mar=c(3,4.5,1.5,3))
        colnames(bottom) <- NULL
        barplot(-bottom, main="", axes=FALSE, legend.text=FALSE, ylab=bottomLab, ...)
        if(is.numeric(cutLineColumn)) abline(v=cutLineColumn*1.2+0.1, lwd=2, lty=2, col="red")
        axis(2)
        grid(col=gridCol, nx=0, ny=NULL)
    }
    
}



###############################################################################################################################################
##### These functions are the three main steps of the pipeline (removing artefacts, estimate the elongation size, generate the 'piled' vectors)
###############################################################################################################################################
#
# Each one is designed to be called on one chromosome at a time, the function using them is taking care of reassembling the chromosomes informations after each call
# Some of them can however be used for several chromosomes at a time, it will list the available chromosomes in the provided alignedDataObject, loop on it and return a list instead of an atomic object
# Be carefull NOT TO use part of chromosomes only... would be tricky to handle the results after... and might be meaningless



# Subfunction that plots the thresholds lines (horizontal and vertical) for the graphs relative to reads occupancy (selection made on threshold)
.plotLineThreshold=function(threshold, percentCumulatedReadsInCounts, line.col)
{
    # Get the x coordinates as numeric
    pileHeight <- as.numeric(names(percentCumulatedReadsInCounts))
    
    # Get the bigger value under threshold
    closestHeight <- sapply(lapply(threshold, function(threshold,values){return(values[values<=threshold])}, pileHeight), max)
    percentCorresponding <- percentCumulatedReadsInCounts[as.character(closestHeight)]
    
    
    #plot vertical line
    segments(x0=threshold,y0=0, x1=threshold, y1=percentCorresponding, lty=2, col=line.col)
    #plot horizontal line
    segments(x0=0,y0=percentCorresponding, x1=threshold, y1=percentCorresponding, lty=2, col=line.col)
    
    text(threshold, percentCorresponding,labels=paste(threshold, " -> ", format(percentCorresponding,digits=4, scientific=FALSE), "%", sep=""),adj=c(-0.1,1.2), col=line.col)        
}

# Subfunction that plots the thresholds lines (horizontal and vertical) for the graphs relative to reads occupancy (selection made on percentage of reads to keep)
#    plotLinePercent=function(percent, percentCumulatedReadsInCounts, line.col)
#    {
#        # Get the x coordinates as numeric
#        pileHeight=as.numeric(names(percentCumulatedReadsInCounts))
#        names(pileHeight)=percentCumulatedReadsInCounts
#        
#        # Get the bigger value under threshold
#        closestPercent=sapply(lapply(percent, function(percent,values){return(values[values>=percent])}, percentCumulatedReadsInCounts),min)
#        pileHeightCorresponding=pileHeight[as.character(closestPercent)]
#        
#        
#        #plot vertical line
#        segments(x0=pileHeightCorresponding,y0=0, x1=pileHeightCorresponding, y1=percent, lty=2, col=line.col)
#        #plot horizontal line
#        segments(x0=0,y0=percent, x1=pileHeightCorresponding, y1=percent, lty=2, col=line.col)
#        
#        text(pileHeightCorresponding, percent,labels=paste(format(percent,digits=4, scientific=FALSE), "%", " -> ", pileHeightCorresponding, " -> ", format(closestPercent,digits=4, scientific=FALSE), "%", sep=""),adj=c(-0.1,1.2), col=line.col)
#        
#    }


# This function will plot two graphs relative to reads distribution along the chromosome (This function doesn't take the strand in account, plase split it before if separate analysis required)
# It will also return as a list the computed parameters that might be necessary for the artefact removal
.plotReadOccupancyStats=function(positionTAB, totalNbReads, thresholds, suffixTitle="")
{
    
    # Counting the number of reads aligned to each position
    positionCOUNT <- table(positionTAB)
    tableCount <- table(positionCOUNT)
    
    
    ### Plotting stats about piles
    
    ## 1 piles height distribution
    
    suppressWarnings(plot(tableCount, log="x", main=paste("Piles height distribution - ",suffixTitle, " - Total nb of reads : ", totalNbReads,sep=""), xlab=paste("Nb of piled reads (pile height) - max=",max(positionCOUNT),sep=""), ylab="Occurences (number of piles)"))
    
    ## 2 Getting stats about nb of reads according to eventual threshold
    # Even if we have these stats, we keep on extracting artefacts (the actual index of reads to remove) for all "stat" threshold in the calling function because we want to know more about paired-ends later
    
    # Was used while plotting two strands on the same graph 
    ## Get the maximum pile size to plot the graph (x axis)
    #maxHeightPile=max(positionCOUNT)
    ## Actual plotting
    #plot(NULL, ylim=c(0.1,100), xlim=c(0.1, maxHeightPile), main="", xlab="",ylab="")
    
    ## Plotting the curves    
    #for(currentStrand in strandValues)
    #{
    # Compute the percentage of reads remaining for each threshold
    nbReadsInCounts <- tableCount*as.numeric(names(tableCount))
    cumulatedReadsInCounts <- cumsum(nbReadsInCounts)
    percentCumulatedReadsInCounts <- (cumulatedReadsInCounts/max(cumulatedReadsInCounts))*100
    
    #points(as.numeric(names(percentCumulatedReadsInCounts)), percentCumulatedReadsInCounts,type="s", col=as.numeric(factor(currentStrand, levels=strandValues)))        
    #}
    
    # Plot the curves
    plot(as.numeric(names(percentCumulatedReadsInCounts)), percentCumulatedReadsInCounts,type="s", main=paste("Proportion of remaining reads by threshold - ", suffixTitle, sep=""), ylab="Remaining reads (%)" , xlab="Threshold (pile height) for artefacts")
    
    # Plot the lines corresponding to selected thresholds
    .plotLineThreshold(thresholds, percentCumulatedReadsInCounts, rainbow(length(thresholds)))
    
    # Return the extracted and computed vectors that will be useful for the artefact removal steps 
    return(NULL)
}


############################################
# "getArtefactsIndexes" will check for the artefacts and return the indexes to remove FROM THE OBJECT USED FOR THE CALL 
# (has to be translated to original indexes if you send a subselection)
# If 'remove=TRUE' the function will return the object used for the call minus the reads concerned
# If 'remove=FALSE' the function will return the index to remove in the object alignedDataObject used for the call
# The second option is better for multithreading, specially for 'multicore' because it avoid copying all the data through the pipe.
# Text messages are concatenated and showed in a unique cat in case of multithreading, it helps not to mix all messages from different chromosomes.
# THE INPUT SHOULD NEVER BE A LIST OF AlignedData OBJECT SPLITTED BY CHROMOSOMES, IT MUST BE AN AlignedData OBJECT ALREADY SPLITTED : ONE CHROMOSOME AT A TIME)
getArtefactsIndexes <- function(alignedDataObject, expName, thresholdToUse=1, thresholdForStats=c(1:5,10,20,50,100), resultFolder=".")
{
    
    # Handle the pass-by-environment case
    if(is.environment(alignedDataObject))
    {
        alignedDataObject <- alignedDataObject$value
    }
    
    
    # Check the chromosomes information in the object
    currentChr <- unname(unique(seqnames(alignedDataObject)))
    
    if(length(currentChr)>1)
    {
        stop("getArtefactIndexes expects an AlignedData object with information on a unique chromosome (seqnames)... Please split your object before calling this function.")
    }
    
    cat(paste(currentChr, "... ", sep=""))
    
    # Basic tests on function arguments to exclude obvious misuses
    if(!(is.numeric(thresholdToUse) || is.numeric(thresholdForStats))) stop("The getArtefactsIndexes function requires at least a threshold to use or a vector of thresholds for stats on artefacts")
    if(length(thresholdToUse)>1) stop("The thresholdToUse argument of function getArtefactsIndexes must be an atomic numeric value, NA, or NULL")
    
    
    
    # Concatenate the thresholds to treat them all at the same time in loops
    thresholds <- as.numeric(unique(na.omit(c(thresholdToUse,thresholdForStats))))
    
    # Will receive all the intermediary or final information to be kept along the analysis or to be returned to calling environment
    listResults <- list()
    
    # Used for the loops
    #strandValues <- c("+","-")
    strandValues <- unname(unique(strand(alignedDataObject)))
    
    
    # Preparing the picture file and layout in which will be plotted the 5 figures (two for statistics on piles for each strand and the last one summarizing removed artefacts)
    skipPlotting <- FALSE
    if(length(strandValues)==2)
    {
        pdf(file=file.path(resultFolder, paste(expName, "_pileStats_",currentChr,".pdf",sep="")), width=10, height=10)
        layout(matrix(c(1,1,3,3,2,4,2,4,5,5,6,6), ncol=2, byrow=TRUE)) # Create a layout for figures organization
    } else if(length(strandValues)==1)
    {
        pdf(file=file.path(resultFolder, paste(expName, "_pileStats_",currentChr,".pdf",sep="")), width=10, height=10)
        layout(matrix(c(1:3), ncol=1)) # Create a layout for figures organization
    } else
    {
        skipPlotting <- TRUE
    }
    
    
    #### 1- This block will select the reads subpopulations (strands) and count their occurences by position (distribution of pile heights)
    for(currentStrand in strandValues)
    {
        indexSelectionFromTotal <- which( (seqnames(alignedDataObject)==currentChr) & (strand(alignedDataObject)==currentStrand) )
        
        # Plot the reads occupancy distribution and store the countings for each strand in the list
        # Getting the coordinates of the reads
        positionTAB <- position(alignedDataObject[indexSelectionFromTotal])
        listResults[[currentStrand]] <- list("positionTAB"=positionTAB, "positionCOUNT"=table(positionTAB))
        if(!skipPlotting)
        {
            .plotReadOccupancyStats(positionTAB, length(alignedDataObject[indexSelectionFromTotal]), thresholds=thresholds, suffixTitle=paste("strand '", currentStrand,"'",sep=""))
        }
        listResults[[currentStrand]][["indexSelectionFromTotal"]] <- indexSelectionFromTotal
    }
    
    # Eventual block of automatic 'selection of' threshold or 'fine tuning' based on the precomputed distributions (see previous block)
    # This is the main reason why the two loops are separated and why there is this mechanism of 'listResult' to keep information between loops
    
    #### 2- This block finds the actual artefacts indexes
    for(currentStrand in strandValues)
    {
        # Will help to keep trace of the threshold for the future lapplys
        names(thresholds) <- thresholds
        
        # Getting the coordinates which have more than 'threshold' reads aligned (and for which we should keep only one copy)
        #positionsRepeatedMoreThanThr=names(listResults[[currentStrand]][["positionCOUNT"]][listResults[[currentStrand]][["positionCOUNT"]]>threshold])
        #This one makes the job on several threshold, and return the result as a list, each element being the result of a threshold
        listResults[[currentStrand]][["positionsRepeatedMoreThanThr"]] <- lapply(lapply(lapply(thresholds,"<",listResults[[currentStrand]][["positionCOUNT"]]), function(x){return(listResults[[currentStrand]][["positionCOUNT"]][x])}),names)
        
        #messages <- paste(messages, "\n    listing unique IDs to remove... nb of positions concerned :", length(listResults[[currentStrand]][["positionsRepeatedMoreThanThr"]]))
        
        
        # Get index in PositionTAB of all reads that are repeated more than 'threshold' times (the ones we gonna REMOVE)
        #indexRepeatedInPositionTAB=(listResults[[currentStrand]][["positionTAB"]] %in% positionsRepeatedMoreThanThr[[as.character(threshold)]])
        #This one makes the job on several threshold, and return the result as a list, each element being the result of a threshold
        listResults[[currentStrand]][["indexRepeatedInPositionTAB"]] <- lapply(listResults[[currentStrand]][["positionsRepeatedMoreThanThr"]],function(positionsRepeatedMoreThanThr){return(listResults[[currentStrand]][["positionTAB"]] %in% positionsRepeatedMoreThanThr)})
        
        
        
        # In case of paired end, we cannot take a decision on which kind of read to keep among the ones for which both are in artefactual regions or not 
        # Indirectly, it avoids to take the risk of breaking the reads/mates symmetry
        if(!pairedEnds(alignedDataObject))
        {
            # Get index in PositionTAB of the the FIRST read matching each repeated position (the only one we gonna KEEP for each repeated position) 
            # indexFirstRepeated <- match(listResults[[currentStrand]][["positionsRepeatedMoreThanThr"]][[as.character(thresholdToUse)]], positionTAB)
            # Do it on all thresholds
            indexFirstRepeated <- lapply(listResults[[currentStrand]][["positionsRepeatedMoreThanThr"]], match, listResults[[currentStrand]][["positionTAB"]])
            
            # We want to keep at least one copy for the reads repeated more than thr, we have to remove all repeats BUT the first occurence for each concerned coordinate
            #indexRepeatedInPositionTAB[indexFirstRepeated] <- FALSE
            # Do it on all thresholds
            listResults[[currentStrand]][["indexRepeatedInPositionTAB"]] <- mapply(function(indexRepeatedInPositionTAB, indexFirstRepeated)
                    {
                        indexRepeatedInPositionTAB[indexFirstRepeated] <- FALSE 
                        return(indexRepeatedInPositionTAB)
                    },
                    listResults[[currentStrand]][["indexRepeatedInPositionTAB"]],
                    indexFirstRepeated, SIMPLIFY=FALSE)
            
        }
        
        #messages <- paste(messages, "\n    list to remove",ifelse(pairedEnds,":", "(keep the first occurence) :"), sum(indexRepeatedInPositionTAB))
        
        # Translating the indexes to the ones in the object having all reads (+ and -)
        listResults[[currentStrand]][["indexesToRemoveFromTotal"]] <- lapply(listResults[[currentStrand]][["indexRepeatedInPositionTAB"]], function(indexRepeatedInPositionTAB){return(listResults[[currentStrand]][["indexSelectionFromTotal"]][indexRepeatedInPositionTAB])})
    }
    
    
    #### 3- Once artefacts were computed on both strands, check mate pairs of reads that were considered as "to be removed"
    if(pairedEnds(alignedDataObject))
    {
        for(currentStrand in strandValues)
        {
            otherStrand <- strandValues[!strandValues %in% currentStrand]
            
            
            ### First in pair (FIP)
            # get the "to remove" reads marked as first in pair
#            firsInPairToRemoveFromTotal <- original_indexesToRemoveFromTotal[which(as.logical(bitAnd(64,flag[original_indexesToRemoveFromTotal])))] # first in pair
#            associatedMate <- firsInPairToRemoveFromTotal+1 # This works because we work on an object SORTED by pairs !!!
            # Do it on all thresholds
            FIPindex_in_indexesToRemoveFromTotal <- lapply(listResults[[currentStrand]][["indexesToRemoveFromTotal"]], function(indexesToRemoveFromTotal) {return(which(as.logical(bitAnd(64,flag(alignedDataObject)[indexesToRemoveFromTotal]))))})
            firstInPairToRemoveFromTotal <- mapply("[", listResults[[currentStrand]][["indexesToRemoveFromTotal"]], FIPindex_in_indexesToRemoveFromTotal, SIMPLIFY=FALSE)
            associatedMate <- lapply(firstInPairToRemoveFromTotal, "+", 1) # This works because we work on an object SORTED by pairs !!!
            
            
            # Get the indexes of reads that have their mate already considered as artefacts (for statistics)
            # We search it ON THE OPPOSITE STRAND because by definition pairs have reads on opposite strands
            alreadyConsideredArtefactFIPmate <- mapply("%in%", associatedMate, listResults[[otherStrand]][["indexesToRemoveFromTotal"]], SIMPLIFY=FALSE)
            
            
            ### Second in pair (SIP)
            # get the "to remove" reads marked as second in pair
#            secondInPairToRemoveFromTotal <- original_indexesToRemoveFromTotal[which(as.logical(bitAnd(128,flag[original_indexesToRemoveFromTotal])))] #second in pair
#            associatedMate <- secondInPairToRemoveFromTotal-1 # This works because we work on an object SORTED by pairs !!!
            # Do it on all thresholds
            SIPindex_in_indexesToRemoveFromTotal <- lapply(listResults[[currentStrand]][["indexesToRemoveFromTotal"]], function(indexesToRemoveFromTotal) {return(which(as.logical(bitAnd(128,flag(alignedDataObject)[indexesToRemoveFromTotal]))))})
            secondInPairToRemoveFromTotal <- mapply("[", listResults[[currentStrand]][["indexesToRemoveFromTotal"]], SIPindex_in_indexesToRemoveFromTotal, SIMPLIFY=FALSE)
            associatedMate <- lapply(secondInPairToRemoveFromTotal, "-", 1) # This works because we work on an object SORTED by pairs !!!
            
            # Get the indexes of reads that have their mate already considered as artefacts (for statistics)
            # We search it ON THE OPPOSITE STRAND because by definition pairs have reads on opposite strands
            alreadyConsideredArtefactSIPmate <- mapply("%in%", associatedMate, listResults[[otherStrand]][["indexesToRemoveFromTotal"]], SIMPLIFY=FALSE)
            
            
            ### Process them
            
            # Merge
            #listResults[[currentStrand]][["indexesToRemoveFromTotal_IsPairedArtefact"]] <- mapply(c,alreadyConsideredArtefactFIPmate,alreadyConsideredArtefactSIPmate)
            # Let's try to get them in the same order as 'indexesToRemoveFromTotal', first create a vector of FALSE and then put as TRUE the ones which are pairs (retranslate indexes back, from FIP/SIP to indexesToRemoveFromTotal using the index_in_indexesToRemoveFromTotal)
            listResults[[currentStrand]][["indexesToRemoveFromTotal_isPairedArtefact"]] <- lapply(listResults[[currentStrand]][["indexesToRemoveFromTotal"]], function(x){return(logical(length(x)))}) # Create a boolean vector all FALSE
            
            # FIP
            listResults[[currentStrand]][["indexesToRemoveFromTotal_isPairedArtefact"]] <- mapply(function(indexesToRemoveFromTotal_isPairedArtefact, alreadyConsideredArtefactFIPmate, FIPindex_in_indexesToRemoveFromTotal)
                    {
                        indexesToRemoveFromTotal_isPairedArtefact[FIPindex_in_indexesToRemoveFromTotal[alreadyConsideredArtefactFIPmate]] <- TRUE
                        return(indexesToRemoveFromTotal_isPairedArtefact)
                    },
                    listResults[[currentStrand]][["indexesToRemoveFromTotal_isPairedArtefact"]], alreadyConsideredArtefactFIPmate,FIPindex_in_indexesToRemoveFromTotal)
            
            # SIP
            listResults[[currentStrand]][["indexesToRemoveFromTotal_isPairedArtefact"]] <- mapply(function(indexesToRemoveFromTotal_isPairedArtefact, alreadyConsideredArtefactSIPmate, SIPindex_in_indexesToRemoveFromTotal)
                    {
                        indexesToRemoveFromTotal_isPairedArtefact[SIPindex_in_indexesToRemoveFromTotal[alreadyConsideredArtefactSIPmate]] <- TRUE
                        return(indexesToRemoveFromTotal_isPairedArtefact)
                    },
                    listResults[[currentStrand]][["indexesToRemoveFromTotal_isPairedArtefact"]], alreadyConsideredArtefactSIPmate, SIPindex_in_indexesToRemoveFromTotal)
            
            # This SHOULD NEVER HAPPEN, or this would reveal an algorithmic failure or a file sorting problem
            #if(!all(ArtefactWithOrphanMate %in% indexesToRemoveFromTotal)) stop("Inconsistency while searching for artefacts with orphan mate (some reads considered as artefacts with orphan mates don't appear in artefacts list before removal), check that your file is properly sorted, first in pair and second in pair should be interlaced")
            
            
            # These SHOULD NEVER BE TRUE, or this would reveal an algorithmic failure or a file sorting problem
            #if(any(orphanMateFromArtefact %in% indexesToRemoveFromTotal)) stop("Inconsistency while searching for orphan mate artefacts (the mate supposed orphan already appears in the artefacts list), check that your file is properly sorted, first in pair and second in pair should be interlaced")
            #if(any(duplicated(c(ArtefactWithOrphanMate,orphanMateFromArtefact,indexesToRemoveFromTotal)))) stop("Inconsistency while searching for artefacts with orphan mate (some artefacts seem to appear in several exclusive lists), check that your file is properly sorted, first in pair and second in pair should be interlaced")
            
        } # /currentStrand (+/-)
    } # /If pairedEnds check mates
    
    
    
    
    ### plotting stats about removed artefacts (by threshold) (5th plot, because the two first are made on each strand)
    if(length(strandValues)==2)
    {
        # Preparing to plot a back to back barplot (one per strand)
        top <- NULL
        bottom <- NULL
        
        # Plot artefacts statistics by thresholds
        if(pairedEnds(alignedDataObject))
        {
            # Top
            artefactToRemoveByThr <- sapply(listResults[["+"]][["indexesToRemoveFromTotal"]],length)
            pairedArtefactsByThr <- sapply(listResults[["+"]][["indexesToRemoveFromTotal_isPairedArtefact"]],sum)
            orphanArtefactToRemoveByThr <- artefactToRemoveByThr-pairedArtefactsByThr
            
            top <- rbind(pairedArtefactsByThr, orphanArtefactToRemoveByThr)
            rownames(top) <- c("Paired artefacts","Artefacts with orphan")
            
            # Bottom
            artefactToRemoveByThr <- sapply(listResults[["-"]][["indexesToRemoveFromTotal"]],length)
            pairedArtefactsByThr <- sapply(listResults[["-"]][["indexesToRemoveFromTotal_isPairedArtefact"]],sum)
            orphanArtefactToRemoveByThr <- artefactToRemoveByThr-pairedArtefactsByThr
            
            bottom <- rbind(pairedArtefactsByThr, orphanArtefactToRemoveByThr)
            rownames(bottom) <- c("Paired artefacts","Artefacts with orphan")
        }
        else
        {
            top <- sapply(listResults[["+"]][["indexesToRemoveFromTotal"]],length)
            bottom <- sapply(listResults[["-"]][["indexesToRemoveFromTotal"]],length)
        }
        
        .plotBackToBackBarplot(top, bottom, topLab="Positive strand (Reads)", bottomLab="Negative strand (Reads)", main="Artefacts distribution", cutLineColumn=1, outsideLayout=TRUE, plotLegend=pairedEnds(alignedDataObject), beside=FALSE)
    }
        
    if(!skipPlotting)
    {
        ### Finish plotting stats
        dev.off()
    }    
    
    ### Select only interesting results to report to global function (Try not to overflow memory especially for 'multithreading' because all the data will be serialized and passed through a pipe, serialization fails for too large objects) 
    
    # Select only results relative to information extracted from strands
    listResults <- listResults[strandValues]
    for(currentStrand in strandValues)
    {
        # select only indexesToRemoveFromTotal and indexesToRemoveFromTotal_isPairedArtefact
        listResults[[currentStrand]] <- listResults[[currentStrand]][c("indexesToRemoveFromTotal","indexesToRemoveFromTotal_isPairedArtefact")]
        # return only results about the thresholdToUse
        listResults[[currentStrand]] <- lapply(listResults[[currentStrand]],"[", as.character(thresholdToUse))
    }
    
    return(listResults)
    
}



############################################
# "estimateElongationSize" will return the most probable elongation size for the provided chromosomes in 'alignedDataObject'
estimateElongationSize <- function(alignedDataObject, expName, stepMin=50, stepMax=450, stepBy=10, averageReadSize=NA, resultFolder=".")
{
    
    # Handle the pass-by-environment case
    if(is.environment(alignedDataObject))
    {
        alignedDataObject <- alignedDataObject$value
    }
    
    
    # Check if there is one or several chromosomes provided in alignedDataObject
    chromosomeNames <- unname(unique(seqnames(alignedDataObject)))
    
    # List of file names created (if there is several chromosomes, it will be a list)
    resultingElongation <- NULL
    if(length(chromosomeNames)>1)
    {
        resultingElongation <- list()
    }
    
    
    for(currentChr in chromosomeNames)
    {
        
        cat(paste(currentChr, "... ", sep=""))
        
        #selection of the reads on strand +
        alignedDataPlus <- alignedDataObject[seqnames(alignedDataObject)==currentChr  &  strand(alignedDataObject)=="+"]
        
        #selection of the reads on strand -
        alignedDataMinus <- alignedDataObject[seqnames(alignedDataObject)==currentChr  &  strand(alignedDataObject)=="-"]
        
        
        # Piling the reads in 'per coordinate' values
        piledP <- .newPileup(position(alignedDataPlus), qwidth(alignedDataPlus))
        
        piledM <- .newPileup(position(alignedDataMinus), qwidth(alignedDataMinus))
        
        piledP <- as.numeric(piledP)
        piledM <- as.numeric(piledM)
        
        stepShift <- seq(stepMin,stepMax,stepBy)
        
        # Computing the overlap score for each shifting step
        res <- elongationEstimation(piledP,piledM,stepShift)
        
        maxScoredShift <- stepShift[which.max(res)]
        
        pdf(file=file.path(resultFolder, paste(expName, "_extensionEstimation_",currentChr,"_",maxScoredShift+(averageReadSize),".pdf",sep="")))
        
        if(is.na(averageReadSize)) # in case a global value has not been provided, get it from the current chromosome
        {
            averageReadSize <- trunc(mean(qwidth(alignedDataObject)))
        }
        
        plot(x=stepShift,y=res,main=paste("Chromosome",currentChr,"size :",maxScoredShift+(averageReadSize),sep=" "))
        abline(v=maxScoredShift)
        
        dev.off()
        
        # Stores the score of the current chromosome in this tab indexed by name of the chromosome
        elongationSizeEstimation <- maxScoredShift+(averageReadSize)
        
        rm(alignedDataPlus, alignedDataMinus, piledP, piledM, res)
        invisible(gc())
        
        
        if(length(chromosomeNames)>1)
        {
            resultingElongation[[currentChr]] <- elongationSizeEstimation
        }
        else
        {
            resultingElongation <- elongationSizeEstimation
        }
    }
    
    return(resultingElongation)
    
}



############################################
# "generatePiled" will generate a vector of scores for Coordinates spanning the chromosome(s) provided in 'alignedDataObject'
# It can handle several situations, basically all possible combination of : Paired-end/Single-end, UniReads/MultiReads (weigthed), MidPoint/Normal piling, Elongation/No Elongation, ReadSize Specified/Not Specified 
generatePiled <- function(alignedDataObject, elongationSize, averageReadSize, midPoint=FALSE)
{
    
    # Handle the pass-by-environment case
    if(is.environment(alignedDataObject))
    {
        alignedDataObject <- alignedDataObject$value
    }
    
    
    if(!((is.numeric(elongationSize) && elongationSize>=0 && length(elongationSize)==1) || is.na(elongationSize))) stop("elongationSize parameter must be an atomic positive number or NA value")
    
    if(!is.na(elongationSize))
    {
        if(pairedEnds(alignedDataObject) && (elongationSize!=0))
        {
            warning("A non-zero shifting (elongation) value has been specified to pileup module but the data is claimed to be paired ends, ignoring parameter (NA expected). Consider declaring your data as single-ended before asking for a custom elongation size.")
            elongationSize <- NA
        }
    }
    else if(!pairedEnds(alignedDataObject))
    {
        stop("A non valid elongation has been specified for a dataset declared as single-ended experiment (NA instead of positive or 0 value expected)")    
    }
    
    # Check if there is one or several chromosomes provided in alignedDataObject
    chromosomeNames <- unname(unique(seqnames(alignedDataObject)))
    
    # List of file names created (if there is several chromosomes, it will be a list)
    resultingPiled <- NULL
    if(length(chromosomeNames)>1)
    {
        resultingPiled <- list()
    }
    
    
    for(currentChr in chromosomeNames)
    {
        alnCurrentChrom <- alignedDataObject[seqnames(alignedDataObject)==currentChr]
        
        cat(paste(currentChr, " (",length(alnCurrentChrom)," reads)","... ", sep=""))
        
        #cat("\n Piling chromosome :", currentChr, "- reads number :",length(alnCurrentChrom), "-", ifelse(is.na(elongationSize),"automatic",paste(elongationSize,"bp")), ifelse(midPoint,"shifting","extension"),"...")
        
        # Removing one read of the pair to avoid double representation of the signal while elongating/shifting
        #if(pairedEnds(alnCurrentChrom) && is.na(elongationSize))
        #{
        #    alnCurrentChrom <- alnCurrentChrom[strand(alnCurrentChrom)=="+"]
        #}
        
        
        ### Set readLengthParam value according to what information is given (Transparently handles multiread by taking a decision on readLength param)
        readLengthParam <- NULL
        
        if(length(qwidth(alnCurrentChrom))>0) # If provided in the object
        {
            readLengthParam <- qwidth(alnCurrentChrom) # Use the real one
        }
        else
        {
            if(!is.numeric(averageReadSize) && averageReadSize>0 && (length(averageReadSize)==1)) 
            {
                stop("No information about read length in object and no valid average read size specified (positive number)")
            }
            readLengthParam <- averageReadSize
        }
        
        
        ### Set fragLengthParam value and do specific treatment according to what information is given
        fragLengthParam <- NULL
        
        
        
        if(midPoint)
        {
            fragLengthParam <- 1
            
            if(pairedEnds(alnCurrentChrom) && is.na(elongationSize)) # Shifting reads by isize/2
            {
                # We only have one read per pair left, the one on the positive strand because we removed the other one to avoid double representation of the signal
                position(alnCurrentChrom) <- as.integer(position(alnCurrentChrom)+trunc(isize(alnCurrentChrom)/2))
                
                #position(alnCurrentChrom[strand(alnCurrentChrom)=="+"])=position(alnCurrentChrom[strand(alnCurrentChrom)=="+"])+trunc(isize(alnCurrentChrom[strand(alnCurrentChrom)=="+"])/2)
                #position(alnCurrentChrom[strand(alnCurrentChrom)=="-"])=position(alnCurrentChrom[strand(alnCurrentChrom)=="-"])-trunc(isize(alnCurrentChrom[strand(alnCurrentChrom)=="-"])/2)
            }
            #else # either (SE & whatever elongation) or (PE & elongation 0) because PE and elongation!=0 ignored at the beginning of the function
            else if((!pairedEnds(alnCurrentChrom)) && (elongationSize!=0)) # Limit it to SE and no need to go here for elongation/shifting 0 (elongationSize can't be NA here because of the tests @ the beginning of the function)
            {
                position(alnCurrentChrom) <- as.integer(position(alnCurrentChrom)+ifelse(strand(alnCurrentChrom)=="+",trunc(elongationSize/2), -trunc(elongationSize/2)))
                #position(alnCurrentChrom[strand(alnCurrentChrom)=="+"])=as.integer(position(alnCurrentChrom[strand(alnCurrentChrom)=="+"])+trunc(elongationSize/2))
                #position(alnCurrentChrom[strand(alnCurrentChrom)=="-"])=as.integer(position(alnCurrentChrom[strand(alnCurrentChrom)=="-"])-trunc(elongationSize/2))
            }
        }
        else
        {
            fragLengthParam <- readLengthParam
            
            if(pairedEnds(alnCurrentChrom) && is.na(elongationSize)) 
            {
                fragLengthParam <- isize(alnCurrentChrom)
            }
            else if((!pairedEnds(alnCurrentChrom)) && (elongationSize!=0))
            {
                fragLengthParam <- elongationSize
            }
        }
        
        
        
        if(is.null(fragLengthParam) || is.null(readLengthParam)) 
        {
            stop("The piling module failed to choose a behaviour from data information, please check that you give all necessary informations (see documentation)")
        }
        
        ### Setting weight parameter as 1 if not specified
        weightParam <- if(length(weight(alnCurrentChrom))>0) weight(alnCurrentChrom) else rep(1.0, length(alnCurrentChrom))
        
        ### PILING
        
        # In case of elongation size specified, it overrides the readLength parameter
        # readLEngth is just used to shift the start coordinate for negative strands
        # so it's better to put the read size when elongationsize=0 but it's still 
        # possible to put less in case needed (nucelosomes start for instance)
        
        
#        piledRLE = .newPileup(    start=position(alnCurrentChrom),
#                                fragLength=fragLengthParam,
#                                dir=strand(alnCurrentChrom),
#                                readLength=qwidth(alnCurrentChrom))
        
        
        piledRLE <- Rle(pileupDouble(start=position(alnCurrentChrom),
                        fragLength=fragLengthParam,
                        dir=strand(alnCurrentChrom),
                        readLength=readLengthParam, 
                        weight=weightParam))
        
        # Set the chromosome information in the Rle object metadata in order to be able to retrieve it in further processing (writing files)    
        metadata(piledRLE) <- list("chr"=currentChr)
        
        # Setting the two vectors at the same size
#        if(length(piled)>length(piledMulti))
#        {
#            length(piledMulti)=length(piled)
#            piledMulti[is.na(piledMulti)]=0
#        }
#        else
#        {
#            length(piled)=length(piledMulti)
#            piled[is.na(piled)]=0
#        }
#
#        # Adding the score of uni- and multi- reads
#        piledRLE=Rle(piled+piledMulti)
        
        
        
        
        if(length(chromosomeNames)>1)
        {
            resultingPiled[[currentChr]] <- piledRLE
        }
        else
        {
            resultingPiled <- piledRLE
        }
    }
    
    return(resultingPiled)
    
}








.writeWIGvs_chr <- function(piledRleData, baseFileName, resultFolder=".", compatibilityOutputWIG=FALSE)
{
    
    currentChr <- metadata(piledRleData)[["chr"]]
    
    cat(paste(currentChr, "... ", sep=""))
    
    # Creating the file and writing the track descriptions
    outputFileName <- file.path(resultFolder, paste("TEMP_WIGvs_", baseFileName, ".", currentChr, sep=""))
    
    # Binary mode to avoid conversion of \n for specific platforms
    fileCon <- file(outputFileName, open="wb")

    # Repeat track line for every chromosome only in compatibility mode (non standard WIG format) 
    if(compatibilityOutputWIG)
    {
        trackLine <- paste("track type=wiggle_0 name=\"", baseFileName, "\"", sep="")
        writeLines(trackLine, con=fileCon, sep="\n")
    }
    
    descLine <- paste("variableStep chrom=", currentChr, sep="")
    writeLines(descLine, con=fileCon, sep="\n")
    
    # Using the RLE to compute the coordinates for variable steps
    
    # Modify the rle to coordinates by summing it and adjust to make it fit (+1)
    startCoord <- diffinv(runLength(piledRleData))+1
    startCoord <- format(startCoord[1:(length(startCoord)-1)], scientific=FALSE, trim=TRUE)
    score <- format(runValue(piledRleData), scientific=FALSE, digits=10, trim=TRUE, drop0trailing = TRUE)
    
    # Writing to file
    writeLines(paste(startCoord, score, sep=" "), con=fileCon, sep="\n")
    
    close(fileCon)
    
    
    resultingFiles <- outputFileName
    
    return(resultingFiles)
    
}


.writeBED_chr <- function(piledRleData, baseFileName, resultFolder=".")
{
    
    currentChr <- metadata(piledRleData)[["chr"]]
    
    cat(paste(currentChr, "... ", sep=""))
    
    # Creating the file and writing the track descriptions
    outputFileName <- file.path(resultFolder, paste("TEMP_BED_", baseFileName, ".", currentChr, sep=""))
    
    # Binary mode to avoid conversion of \n for specific platforms
    fileCon <- file(outputFileName, open="wb")

    
    # Using the RLE to compute the coordinates for variable steps
    
    # Modify the rle to coordinates by summing it (BEDgraph coordinates are O-based as opposed to wiggle ones)
    coordBase <- diffinv(runLength(piledRleData))
    endCoord <- format(coordBase[2:(length(coordBase))], scientific=FALSE, trim=TRUE)
    startCoord <- format(coordBase[1:(length(coordBase)-1)], scientific=FALSE, trim=TRUE)
    score <- format(runValue(piledRleData), scientific=FALSE, digits=10, trim=TRUE, drop0trailing = TRUE)
    
    # Writing to file
    writeLines(paste(currentChr, startCoord, endCoord, score, sep=" "), con=fileCon, sep="\n")
    
    close(fileCon)
    
    
    resultingFiles <- outputFileName
    
    return(resultingFiles)
    
}



.writeWIGfs_chr <- function(binnedDataRle, baseFileName, binSize, resultFolder=".", compatibilityOutputWIG=FALSE)
{
    
    currentChr <- metadata(binnedDataRle)[["chr"]]
    
    cat(paste(currentChr, "... ", sep=""))
    
    # Formatting numbers to remove trailing 0s and avoid scientific notation
    formattedScores <- format(as.numeric(binnedDataRle), scientific=FALSE, digits=10, trim=TRUE, drop0trailing = TRUE)
    
    # Creating the file and writing the track descriptions
    outputFileName <- file.path(resultFolder, paste("TEMP_WIGfs_", baseFileName, ".", currentChr, sep=""))
    
    # Binary mode to avoid conversion of \n for specific platforms
    fileCon <- file(outputFileName, open <- "wb")
    
    # Repeat track line for every chromosome only in compatibility mode (non standard WIG format) 
    if(compatibilityOutputWIG)
    {
        trackLine <- paste("track type=wiggle_0 name=\"", baseFileName, "\"", sep="")
        writeLines(trackLine, con=fileCon, sep="\n")
    }
    
    descLine <- paste("fixedStep chrom=", currentChr, " start=", trunc(binSize/2), " step=", binSize, sep="")
    writeLines(descLine, con=fileCon, sep="\n")
    
    # Writing the binned values
    writeLines(formattedScores, con=fileCon, sep="\n")
    
    close(fileCon)
    
    
    resultingFiles <- outputFileName
    
    return(resultingFiles)
    
}


.writeGFF_chr <- function(binnedDataRle, baseFileName, binSize, resultFolder=".")
{
    currentChr <- metadata(binnedDataRle)[["chr"]]
    
    cat(paste(currentChr, "... ", sep=""))
    
    # Formatting numbers to remove trailing 0s and avoid scientific notation
    formattedScores <- format(as.numeric(binnedDataRle), scientific=FALSE, digits=10, trim=TRUE, drop0trailing = TRUE)
    
    # Creating the file and writing the track descriptions
    outputFileName <- file.path(resultFolder, paste("TEMP_GFF_", baseFileName, ".", currentChr, sep=""))
    
    # Binary mode to avoid conversion of \n for specific platforms
    fileCon <- file(outputFileName, open="wb")
    
    coord <- seq(0, by=binSize, length.out=length(formattedScores))
    
    # Writing the GFF lines (alternative with blocks spanning all the chr)
    #writeLines(paste(currentChr, expName, ".", format(coord, scientific=FALSE, trim=TRUE), format(coord+binSize-1, scientific=FALSE, trim=TRUE), piledBinnedFormatted, ".", ".", format(1:length(piledBinnedFormatted), scientific=FALSE, trim=TRUE), sep="\t"), con=fileCon, sep="\n")
    
    # Writing the GFF lines
    writeLines(paste(currentChr, baseFileName, ".", format(coord+trunc(binSize/2), scientific=FALSE, trim=TRUE), format(coord+trunc(binSize/2), scientific=FALSE, trim=TRUE), formattedScores, ".", ".", paste(currentChr,"_",format(1:length(formattedScores), scientific=FALSE, trim=TRUE), sep=""), sep="\t"), con=fileCon, sep="\n")
    
    close(fileCon)
    
    
    resultingFiles <- outputFileName
    
    return(resultingFiles)
    
}




############################################
# ".readTextFile" will be used to extend the number of file types that readAligned can handle
# It will be used by the function 'readAlignedFiles' above to be able to read GFF and BED files
# It returns an 'AlignedRead' (ShortRead) object with as much information as could be found in the file/arguments
# Default arguments defined for basic BED files.
# If there is no SEQ colum (NA) it will generate a 'fake' sequence filled with Ns based on the size ('end-begin' if no fragmentLength provided)
# columnsIndexes=c(chr=1, begin=2, end=6, strand=5, seq=NA, NULL=4, NULL=3)
.readTextFile <- function(fileName, chrToSelect=NA, columnsIndexes=c(chr=1, begin=2, end=3, strand=4, seq=NA), separator=c(" ", "\t"), positiveStrand=c("F", "1",  "+"), negativeStrand=c("R", "-1", "-"), readLength=NA, chunckSize=(2^19))
{
    
    if(is.na(readLength) && (sum(is.na(columnsIndexes[c("begin", "end")]))>=1) && is.na(columnsIndexes["seq"])) warning("Input : There is no sequence column in the input file and not enough information to get the reads sizes (begin and end column or 'readLength' argument) --> your reads will be assumed to have length 0... very unlikely to be true...")
    
    seq <- list()
    chr <- list()
    position <- list()
    strand <- list()
    
    # Will specify to scan the type of the data columns
    # Must be in the order specified in columnsIndexes (names MUST match !!!)
    myWhat <- list(chr=character(), begin=integer(), end=integer(), strand=character(), seq=character(), NULL=NULL)[names(sort(columnsIndexes))]
    
    con <- file(fileName, open="rt")
    
    # Reading the file by chunks
    while(length((currentRead <- scan(file=con, nmax=chunckSize, what=myWhat)))>0)
    {
#cat(".")
        
        print("1")
        
        chrInSelection <- rep(TRUE,chunckSize)
        
        print("2")
        # Remove the prefix and suffix of the chromosomes
        #currentRead[["chr"]]=gsub(paste(chrPrefix, chrSuffix, sep="|"), "", currentRead[["chr"]]) # Done later now
        print("3")
        if(!is.na(chrToSelect))
        {
            # Select the reads on the desired chromosome
            chrInSelection <- (currentRead[["chr"]] %in% chrToSelect)
        }
        
        if(any(chrInSelection))
        {
            print("4")
            chr <- c(chr,currentRead[["chr"]])
            print("5")
            # Begin and End columns are separated so we can be sure to take the minimum (left) as position of the read (as assumed by readAligned)
            # Even if the strand orientation affect in which wolumn the begin and end are
            beginColumn  <-  as.integer(currentRead[["begin"]])
            endColumn    <-  as.integer(currentRead[["end"]])
            print("6")
            # This will work even if there is no begin or no end column (na.rm)
            leftCoord   <-  mapply(min, beginColumn, endColumn, na.rm=TRUE)
            rightCoord  <-  mapply(max, beginColumn, endColumn, na.rm=TRUE)
            print("7")
            position <- c(position, leftCoord)
            print("8")
            strand <- c(strand,currentRead[["strand"]])
            print("9")
            # If there was no column for the sequence (NA), it will fill it with Ns (nb based on readLength or coordinates)
            if(is.na(columnsIndexes["seq"]))
            {
                seqLengthTAB <- if(is.na(readLength)) (rightCoord-leftCoord) else rep(readLength, sum(chrInSelection))
                print("10")
                # Make a list of unique "N" character and repeat them using seqLengthTAB in mapply, then collapse all Ns for each element using sapply
                replacementSeqs <- DNAStringSet(mapply(function(x,seqLength){paste(rep(x,seqLength),collapse="")},"N", seqLengthTAB))
                print("11")
                seq <- c(seq, replacementSeqs)
                print("12")
            }
            else
            {
                seq <- c(seq,currentRead[["seq"]])
            }
        }
        
    } # /WHILE READLINES
    
    close(con)
    
    seq <- do.call(c, seq)
    chr <- do.call(c, chr)
    position <- do.call(c, position)
    strand <- do.call(c, strand)
    
    # Converting various strand names that can be encountered to the standard factor with appropriate levels
    strand[strand %in% positiveStrand] <- "+"
    strand[strand %in% negativeStrand] <- "-"
    strand <- factor(strand, levels=levels(strand()))
    
    return(ShortRead::AlignedRead(sread=seq, seqnames=as.factor(paste("chr", chr, sep="")), position=position, strand=strand))
}

#system.time((a=.readTextFile()))

############################################
# SAVED COPY
.readTextFile <- function(fileName, chrPrefix="chr", chrSuffix="", chrToSelect=NA, columnsIndexes=c("chr"=1, "begin"=2, "end"=3, "strand"=4, "seq"=NA), separator=c(" ", "\t"), positiveStrand=c("F", "1",  "+"), negativeStrand=c("R", "-1", "-"), readLength=NA, chunckSize=(2^16))
{
    
    if(is.na(readLength) && (sum(is.na(columnsIndexes[c("begin", "end")]))>=1) && is.na(columnsIndexes["seq"])) warning("Input : There is no sequence column in the input file and not enough information to get the reads sizes (begin and end column or 'readLength' argument) --> your reads will be assumed to have length 0... very unlikely to be true...")
    
    seq <- list()
    chr <- list()
    position <- list()
    strand <- list()
    
    con <- file(fileName, open="rt")
    
    # Reading the file by chunks
    while(length((currentLine <- readLines(con=con, n=chunckSize)))>0)
    {
#cat(".")
        chrInSelection <- rep(TRUE,chunckSize)
        
        # Split the line by separators and get them as data.frame or matrix
        splittedLine <- do.call(rbind,strsplit(currentLine, paste("[",paste(separator, collapse="|"), "]", sep="")))
        
        # Remove the prefix and suffix of the chromosomes
        splittedLine[,columnsIndexes["chr"]] <- gsub(paste(chrPrefix, chrSuffix, sep="|"), "", splittedLine[,columnsIndexes["chr"]])
        
        if(!is.na(chrToSelect))
        {
            # Select the reads on the desired chromosome
            chrInSelection <- (splittedLine[,columnsIndexes["chr"]] %in% chrToSelect)
        }
        
        if(any(chrInSelection))
        {
            chr <- c(chr,splittedLine[chrInSelection,columnsIndexes["chr"]])
            
            # Begin and End columns are separated so we can be sure to take the minimum (left) as position of the read (as assumed by readAligned)
            # Even if the strand orientation affect in which wolumn the begin and end are
            beginColumn  <-  as.integer(splittedLine[chrInSelection,columnsIndexes["begin"]])
            endColumn    <-  as.integer(splittedLine[chrInSelection,columnsIndexes["end"]])
            
            # This will work even if there is no begin or no end column (na.rm)
            leftCoord   <-  mapply(min, beginColumn, endColumn, na.rm=TRUE)
            rightCoord  <-  mapply(max, beginColumn, endColumn, na.rm=TRUE)
            
            position <- c(position, leftCoord)
            
            strand <- c(strand,splittedLine[chrInSelection,columnsIndexes["strand"]])
            
            # If there was no column for the sequence (NA), it will fill it with Ns (nb based on readLength or coordinates)
            if(is.na(columnsIndexes["seq"]))
            {
                seqLengthTAB <- if(is.na(readLength)) (rightCoord-leftCoord) else rep(readLength, sum(chrInSelection))
                # Make a list of unique "N" character and repeat them using seqLengthTAB in mapply, then collapse all Ns for each element using sapply
                replacementSeqs <- DNAStringSet(mapply(function(x,seqLength){paste(rep(x,seqLength),collapse="")},"N", seqLengthTAB))
                seq <- c(seq, replacementSeqs)
            }
            else
            {
                seq <- c(seq,splittedLine[chrInSelection,columnsIndexes["seq"]])
            }
        }
        
    } # /WHILE READLINES
    
    close(con)
    
    seq <- do.call(c, seq)
    chr <- do.call(c, chr)
    position <- do.call(c, position)
    strand <- do.call(c, strand)
    
    # Converting various strand names that can be encountered to the standard factor with appropriate levels
    strand[strand %in% positiveStrand] <- "+"
    strand[strand %in% negativeStrand] <- "-"
    
    strand <- factor(strand, levels=levels(strand()))
    
    return(ShortRead::AlignedRead(sread=seq, seqnames=as.factor(paste("chr", chr, sep="")), position=position, strand=strand))
}




############################################################################################################################################
##### This is the main function of the pipeline, calling all the other ones
############################################################################################################################################
#
processPipeline <- function(
        # I/O GENERAL PARAMETERS
        INPUTFilesList,
        resultSubFolder               = "Results_Pasha",
        reportFilesSubFolder          = ifelse(length(resultSubFolder)>1,resultSubFolder[2], "ReportFiles"),
        WIGfs                         = TRUE,
        WIGvs                         = FALSE,
        GFF                           = FALSE,
        BED                           = FALSE,
		BIGWIG						  = FALSE,
        compatibilityOutputWIG        = FALSE,
        # COMPLEX PARAMETERS (ATOMIC OR VECTORS OR LIST OF IT)
        incrArtefactThrEvery          = 7000000,
        binSize                       = 50,
        elongationSize                = NA,
        rangeSelection                = IRanges(0,-1), 
        annotationFilesGFF            = NA, # GFF files
        annotationGenomeFiles         = NA, # path to file or "mm9", "hg19"... 
        # SINGLE PARAMETERS
        elongationEstimationRange     = c(mini=150, maxi=400, by=10),
        rehabilitationStep            = c("orphans","orphansFromArtefacts"),
        removeChrNamesContaining      = "random|hap",
        ignoreInsertsOver             = 500,
        nbCPUs                        = 1,
        keepTemp                      = TRUE, # Keep the intermediary files that led to the final ones (rehab and multi)
        logTofile                     = "./log.txt",
        eraseLog                      = FALSE,
        # LIST PARAMETERS (one element per expName)
        multiLocFilesList             = list()) # A list with experiments names and associated filenames to treat
{
    
    # Checking main arguments
    if(!all(c(is.logical(WIGfs), is.logical(WIGvs), is.logical(GFF), is.logical(BED), is.logical(BIGWIG), is.logical(compatibilityOutputWIG)))) stop("Arguments 'WIGfs', 'WIGvs', 'GFF', 'BED', 'BIGWIG', and 'compatibilityOutputWIG' must be logical (and not NA)...")
    
    if(any(c(is.na(WIGfs),is.na(WIGvs),is.na(GFF),is.na(BED), is.na(BIGWIG), is.na(compatibilityOutputWIG)))) stop("Arguments 'WIGfs', 'WIGvs', 'GFF', 'BED', 'BIGWIG', and 'compatibilityOutputWIG' cannot be NA...")
    
    if(compatibilityOutputWIG)
    {
        warningMessage <- "Argument 'compatibilityOutputWIG' is used for compatibility with previous versions only and generates non-standard WIG files (repeated track line), please consider changing this argument to FALSE..."
        cat("\nWARNING :", warningMessage)
        warning(warningMessage)
    }
    
    if(!is.character(rehabilitationStep)) stop("Argument rehabilitationStep must be a character vector (empty or a vector of values in 'orphans' and 'orphansFromArtefacts')")
    if(length(rehabilitationStep)>0)
    {
        if(!all(rehabilitationStep %in% c("orphans","orphansFromArtefacts"))) stop("Argument rehabilitationStep must be a character vector (empty or a vector of values in 'orphans' and 'orphansFromArtefacts')")    
    }
    
    if(!is.na(ignoreInsertsOver))
    {
        if(!(is.numeric(ignoreInsertsOver) && (length(ignoreInsertsOver)==1) && (ignoreInsertsOver>0))) stop("Argument ignoreInsertsOver must be a strictly positive number or NA")
    }
    
    if(!(is.numeric(elongationEstimationRange) && (length(elongationEstimationRange)==3) && all(names(elongationEstimationRange) %in% c("mini","maxi","by")) && all(elongationEstimationRange>0))) stop("Argument elongationEstimationRange must be a named vector of 3 strictly positive numeric values, named 'mini', 'maxi', 'by'")
    
    if(!(is.character(removeChrNamesContaining) && (length(removeChrNamesContaining)==1))) stop("Argument removeChrNamesContaining must be an atomic (length==1) character string (empty or not)")
    
    if(!(is.logical(keepTemp) && (length(keepTemp)==1))) stop("Argument keepTemp must be an atomic logical value")
    
    if(!(is.logical(eraseLog) && (length(eraseLog)==1))) stop("Argument eraseLog must be an atomic logical value")
    
    if(is.numeric(nbCPUs) && (length(nbCPUs)==1))
    {
        options("mc.cores"=nbCPUs)
        # Taking advantage of setting options for parallel to encapsulate mclapply and mcmapply in a wrapper that automatically desactivate the prescheduling
        mclapply <- function(X, FUN,...) {parallel::mclapply(X, FUN,..., mc.preschedule=FALSE)}
        mcmapply <- function(FUN,...) {parallel::mcmapply(FUN,..., mc.preschedule=FALSE)}
    } else
    {
        stop("nbCPUs argument must be a strictly positive number")
    }
    
    # Checking files description and arguments
    
    inputMandatoryInfo <- c("expName", "fileName", "fileType")
    inputOptionalInfo <- c("chrPrefix", "chrSuffix", "pairedEnds", "midPoint")
    inputOptionalInfo_Default <- c("chrPrefix"="", "chrSuffix"="", "pairedEnds"=FALSE, "midPoint"=FALSE)
    
    # Reading files to treat from a list file
    if(is.character(INPUTFilesList))
    {
        # in case the user provided a file name in which experiments parameters are described
        if(file.exists(INPUTFilesList))
        {
            # Read the file and format its info as a list for processing
            expListFromFile <- read.table(INPUTFilesList, header=TRUE, quote="", sep="\t", stringsAsFactors=FALSE)
            
            # Check that minimal information is present
            if(!all(inputMandatoryInfo %in% colnames(expListFromFile)))
            {
                stop(paste("The input file listing the experiments lacks some mandatory information : ",paste(inputMandatoryInfo, collapse=", "),sep=""))        
            }
            
            # Complete missing columns with default values
            for(currentOption in inputOptionalInfo)
            {
                if(!currentOption %in% colnames(expListFromFile))
                {
                    expListFromFile <- cbind(expListFromFile, currentOption=inputOptionalInfo_Default[currentOption])
                }
            }
            
            # Create a list with filename (used to mapply) and named by experiment names
            INPUTFilesList <- as.list(expListFromFile[["fileName"]])
            names(INPUTFilesList) <- expListFromFile[["expName"]]
            
            # Create the final list with all checked arguments (keeping the experiments names)
            INPUTFilesList <- mapply(function(fileName, fileType, chrPrefix, chrSuffix, pairedEnds, midPoint)
                    {
                        return(list(folderName=dirname(fileName), fileName=basename(fileName), fileType=fileType, chrPrefix=chrPrefix , chrSuffix=chrSuffix, pairedEnds=pairedEnds, midPoint=midPoint))
                    }, 
                    INPUTFilesList, 
                    expListFromFile[["fileType"]], 
                    expListFromFile[["chrPrefix"]], 
                    expListFromFile[["chrSuffix"]], 
                    expListFromFile[["pairedEnds"]], 
                    expListFromFile[["midPoint"]], 
                    SIMPLIFY=FALSE)
        } else
        {
            stop("The file describing experiments does not exist !")
        }
    }
    
    
    # check that the list (either read from file or passed directly as argument) has all mandatory fields
    if(!all(sapply(INPUTFilesList, function(x){c("folderName", "fileName", "fileType", "chrPrefix", "chrSuffix", "pairedEnds", "midPoint") %in% names(x)})))
    {
        stop("Your input file list does not contain all required names. It should contain the following named elements:'folderName', 'fileName', 'fileType', 'chrPrefix', 'chrSuffix', 'pairedEnds', 'midPoint'")
    }
    
    
    
    # Checking that experiment names are all different
    if(any(duplicated(names(INPUTFilesList)))) stop("All experiment names must be different in order to avoid confusing returned result. In case one really wants it, the function should be called twice...")
    
    # Checking if files to process exist
    if(!all(sapply(INPUTFilesList, function(x){file.exists(file.path(x$folderName, x$fileName))}))) stop("At least one filename doesn't refer to an existing file to process...")
    # Same but adapted to regular expression pattern (accept gziped files for instance) this is better since the readAligned function does use regular expression pattern
    # -- Lately te regular expression from readaligned are ignored thanks to "^" and "$" added in the pattern, it avoid to load files with approaching names without knowing it...
    #if(!all(sapply(INPUTFilesList, function(x){length(dir(x$folderName, pattern=x$fileName))==1}))) stop("At least one filename pattern doens't refer to an existing file to process or refers to several files")
    
    # Checking that chrPrefix and chrSuffix are characters tabs of length 1
    if(!all(sapply(INPUTFilesList, function(x){return(is.character(x$chrPrefix) && length(x$chrPrefix)==1)}))) stop("At least one experiment has not chrPrefix defined as atomic character value...")
    if(!all(sapply(INPUTFilesList, function(x){return(is.character(x$chrSuffix) && length(x$chrSuffix)==1)}))) stop("At least one experiment has not chrSuffix defined as atomic character value...")
    
    # Checking that PairedEnds, midPoint and dontSort are boolean tabs of length 1
    if(!all(sapply(INPUTFilesList, function(x){return(is.logical(x$pairedEnds) && length(x$pairedEnds)==1)}))) stop("At least one experiment has invalid 'pairedEnds' value (atomic logical expected)...")
    if(!all(sapply(INPUTFilesList, function(x){return(is.logical(x$midPoint) && length(x$midPoint)==1)}))) stop("At least one experiment has invalid 'midPoint' value (atomic logical expected)...")
    #if(!all(sapply(INPUTFilesList, function(x){return(is.logical(x$dontSort) && length(x$dontSort)==1)}))) stop("At least one experiment has invalid 'dontSort' value (atomic logical expected)...")
    
    if(is.list(multiLocFilesList) && (length(is.list(multiLocFilesList))>0))
    {
        # Checking if files for multiLoc also exist
        if(!all(sapply(multiLocFilesList, file.exists))) stop("At least one of the multiLoc file specified does not exist or is not readable...")
        
        # Checking if files for multiLoc actually refer to an existing experiment !
        if(!all(names(multiLocFilesList) %in% names(INPUTFilesList))) stop("At least one of the multiLoc files list does refer to an experiment name that does not exist in the experiments files list...")
    }
    
    
    #### CHECKING AND REFORMATING IMPORTANT OR COMPLEX PARAMETERS
    # Complex parameters are the ones which can be defined in several ways
    # These ones can be either an atomic value, a vector of values, or a list of previous choices sharing the same names as INPUTFilesList (one per experiment)
    # The '.checkComplexParam' function transforms atomic values or vectors to a final list repeating these parameters for each expriment
    # It will also check for structure (experiment names) and values consistency (validity function) 
    
    incrArtefactThrEvery <- .checkComplexParam(incrArtefactThrEvery, "incrArtefactThrEvery", INPUTFilesList) # NA for no artefact removal, otherwise >=0 for estimation of the threshold based on the number of reads in the experiment, negative value for direct specification of the threshold...
    
    binSize <- .checkComplexParam(binSize, "binSize", INPUTFilesList, validValuesFct=function(x){return( (is.numeric(x)) && (length(x)>0) && all(x>0) )})
    
    elongationSize <- .checkComplexParam(elongationSize, "elongationSize", INPUTFilesList, validValuesFct=function(x){return( ((is.numeric(x) && all(x>=0)) || is.na(x)) && (length(x)>0)) }) # NA for automatic estimation or 0 for no elongation, no negative
    
    rangeSelection <- .checkComplexParam(rangeSelection, "rangeSelection", INPUTFilesList, validValuesFct=function(x){return( ((class(x)=="IRanges") && (length(x)>0) && all(start(x)>=0)) || (is.numeric(x) && (length(x)==1) && (x>1)) )})
    
    # Two validation functions that are a bit more complex are defined outside the call
    .validateAnnotationsFileNames <- function(x)
    {
        # check format
        if(length(x)<1) return(FALSE)
        
        if(length(x)==1 && is.na(x)) return(TRUE)
        # Checking if it's an existing file
        if(is.character(x) && all(file.exists(x))) return(TRUE)
        # Not valid
        return(FALSE)
    }
    
    annotationFilesGFF <- .checkComplexParam(annotationFilesGFF, "annotationFilesGFF", INPUTFilesList, validValuesFct=.validateAnnotationsFileNames)
    
    
    .validateGenomeFileNames <- function(x)
    {
        # check format
        if(length(x)!=1) return(FALSE)
        
        if(is.na(x)) return(TRUE)
        # Checking if it's an existing file
        if(is.character(x) && file.exists(x)) return(TRUE)
        # Otherwise checking if it's a reference to a precomputedReferenceFile in resources ("mm9" ...)
        precomputedReferenceFilesFolder="resources"
        if(sum(nchar(system.file(precomputedReferenceFilesFolder, paste(x, ".ref", sep="") ,package="Pasha")))>0) return(TRUE)
        # Not valid
        return(FALSE)
        
    }
    
    annotationGenomeFiles <- .checkComplexParam(annotationGenomeFiles, "annotationGenomeFiles", INPUTFilesList, validValuesFct=.validateGenomeFileNames)
    

    #pairedEnds=.checkComplexParam(pairedEnds, "pairedEnds", INPUTFilesList, validValuesFct=function(x){return( (is.logical(x)) && (length(x)==1) )})
    #dontSort=.checkComplexParam(dontSort, "dontSort", INPUTFilesList, validValuesFct=function(x){return( (is.logical(x)) && (length(x)==1) )})
    #midPoint=.checkComplexParam(midPoint, "midPoint", INPUTFilesList, validValuesFct=function(x){return( (is.logical(x)) && (length(x)==1) )})
    
    
    ####
    
#
#    if((!is.null(previousReturn)) & (!is.na(previousReturn)))
#    {
#        if(is.character(previousReturn))
#        {
#            trycatch(
#            {
#                load(previousReturn)
#            },
#            error=function(x)
#            {
#                stop("The 'previousReturn' parameter is not a valid R file")
#            })
#            previousReturn=returnList
#        }
#        else if(is.list(previousReturn))
#        {
#            # Has been removed because when a computation crash it might be useful to use it, even if there is not all experiments informations
#            #if(!all(names(INPUTFilesList) %in% names(previousReturn))) stop("The 'previousReturn' argument doens't contain information on all experiments")
#        }
#        else
#        {
#            stop("The 'previousReturn' parameter must be a list or a valid filename")
#        }
#    }
    
    ####
    
    # Flag indicating whether a generallog file has finally been created or not (to release it at the end)
    generalLogFile <- FALSE
    
    if(is.null(logTofile))
    {
        cat("\n\nNo general log file specified...")
    } else if((length(logTofile)==1) & is.character(logTofile) & (nchar(logTofile)>0))# Redirecting all the console messages to the specified log file
    {
        if(file.exists(logTofile) && (!eraseLog)) stop("The general log file (for all experiments) specified already exists, erase or rename it first !")
        sink(file=logTofile, split=TRUE)
        generalLogFile <- TRUE
        cat("\n")
        cat("\nLog to file :", logTofile)
    } else
    {
        stop("Invalid log file parameter (expected NULL or a valid filename)")
    }
    
    #cat("\n")
    #cat("\nCALL :", deparse(match.call))
    
    cat("\n\n")
    
    cat("\nNumber of experiment(s) :", length(INPUTFilesList))
    
    
    ## COMPLEX PARAMETERS
    
    cat("\nArtefacts threshold :") # , if(is.na(incrArtefactThrEvery) | incrArtefactThrEvery<=0) "NotRemoving" else incrArtefactThrEvery)
    mapply(function(expName, incrArtefactThrEverys)
            {
                cat("\n   ",expName,":", replace(incrArtefactThrEverys, is.na(incrArtefactThrEverys),"NotRemoving"))
            }, names(incrArtefactThrEvery), incrArtefactThrEvery)
    
    cat("\nBins size :")# , if(is.na(binSize) | binSize<0) "Not piling" else binSize)
    mapply(function(expName, binSizes)
            {
                cat("\n   ",expName,":", replace(binSizes, is.na(binSizes),"NoPiling"))
            }, names(binSize), binSize)
    
    cat("\nElongation size (shifting if midpoint) :")# , if(is.na(elongationSize) | elongationSize<0) "Automatic" else elongationSize)
    mapply(function(expName, elongationSizes)
            {
                cat("\n   ",expName,":", replace(elongationSizes, is.na(elongationSizes),"Automatic"))
            }, names(elongationSize), elongationSize)
    
    cat("\nSelection of inserts/reads specific range(s) :")# , if(is.na(elongationSize) | elongationSize<0) "Automatic" else elongationSize)
    resNULL <- mapply(function(expName, rangeSelection)
            {
                cat("\n   ",expName,":", if(is.numeric(rangeSelection)) paste(rangeSelection, "equal groups and a group with all") else paste(ifelse(width(rangeSelection)==0,"AllReads",paste(start(rangeSelection),"-", end(rangeSelection),sep="")), collapse=" "))
            }, names(rangeSelection), rangeSelection)
    
    ## EXPERIMENT PARAMETERS
    
    cat("\nSequencing and aligning strategy :")
    mapply(function(expName, expParams)
            {
                cat("\n   ",expName,"declared as :", ifelse(expParams$pairedEnds, "Paired ends","Single end"))
            }, names(INPUTFilesList), INPUTFilesList)
    
    cat("\nElongation and piling strategy :")
    mapply(function(expName, expParams)
            {
                cat("\n   ",expName,"will be piled with", ifelse(expParams$midPoint, "MIDPOINT", "classic"), "strategy")
            }, names(INPUTFilesList), INPUTFilesList)
    
    
#
#cat("\nMultiple location reads file :")
#mapply(function(expName, multiLocFile)
#        {
#            cat("\n   ",expName,":", replace(multiLocFile, is.na(multiLocFile),"No multiloc file"))
#        }, names(multiLocFilesList), multiLocFilesList)
    
    cat("\n\nOutput piled-up format(s)", paste(if(WIGfs) "WIG fixed steps", if(WIGvs) "WIG variable steps", if(GFF) "GFF fixed steps", if(BED) "Bedgraph variable step", sep= " - "), sep=" : ")
    cat("\nNumber of CPUs or Cores to use (if available) :", nbCPUs)
    cat(if(keepTemp) "\nThe temporary output files for each pileup category will be kept" else "\nThe temporary output files for each pileup category will be discarded")
    cat("\nWhen applicable, elongation will be estimated from ", elongationEstimationRange["mini"], "bp to ", elongationEstimationRange["maxi"], "bp, every ", elongationEstimationRange["by"], "bp", sep="")
    
    if(length(rehabilitationStep)>0)
    {
        cat("\nRehabilitation module will be loaded with (if applicable) :", paste(rehabilitationStep, collapse=" - "))
    } else
    {
        cat("\nNo rehabilitation for orphan reads")
    }
    
    if(nchar(removeChrNamesContaining)>0)
    {
        cat("\nChromosomes filter will remove chromosome names (seqnames) containing :", removeChrNamesContaining)
    } else
    {
        cat("\nNot filtering chromosome by names")
    }
    
    if(!is.na(ignoreInsertsOver))
    {
        cat("\nIf applicable paired-end inserts over ",ignoreInsertsOver, "bp will be ignored", sep="")
    } else
    {
        cat("\nNot filtering inserts outliers by size")
    }
    
    
    # Will keep the elements to return for each experiment
    returnList <- list()
    
    for(expName in names(INPUTFilesList))
    {
        
        returnList[[expName]] <- list()
        # List storing the time of execution for the main steps (useful for benchmarking different nb of processors)
        returnList[[expName]][["execTime"]] <- list()
        
        # This name is used to store the progress accross the program loops in order to report statistics for each steps
        programProgressName <- expName
        
        
        # Creating a general folder to store results
        resultFolder <- file.path(INPUTFilesList[[expName]]$folderName, paste(resultSubFolder[1], ifelse(INPUTFilesList[[expName]]$pairedEnds,"_PE","_SE" ), sep=""))
        
        .safeCreateFolder(resultFolder)
        
        
        # Create a timestamp for the local log file
        currentTime <- as.POSIXlt(Sys.time())
        stampPrefix <- format(currentTime, format="%Y_%m_%d_%Hh%M")
        
        cat("\n\nProgram version :",as.character(packageVersion("Pasha")))
        
        
        # Redirecting the log for this experiment specifically in the log folder of the experiment
        localLogToFile <- file.path(resultFolder, paste(expName,"_",stampPrefix,"_Pasha.log",sep=""))
        sink(file=localLogToFile, split=TRUE)
        
        
        cat("\n\nLocal log to file :", localLogToFile)
        
        
        cat("\n\n")
        ts <- paste("##------", date(), "------##")
        cat("\n",ts,"\n")
        title <- paste("------ EXP : ",expName," ------",sep="")
        cat("\n")
        cat(paste(rep("-",nchar(title)),collapse=""),"\n")
        cat(title,"\n")
        cat(paste(rep("-",nchar(title)),collapse=""),"\n")
        
        
        #### READING ALIGNED FILE
        
        cat("\nREADING ALIGNED FILE")
        
        alignedDataObject <- NULL
        
        # If the 'previousReturn' return parameter is specified, try to get the binary aligned object from it (directly from the object or on disk), otherwise reads the original file name
#        if("binaryObject" %in% names(previousReturn[[expName]]))
#        {
#            cat("\n From previous analysis (parameter)")
#            alignedDataList=previousReturn[[expName]][["binaryObject"]]
#        }
#        else if("binaryObjectFileName" %in% names(previousReturn[[expName]]))
#        {
#            cat("\n From previous analysis (file)")
#            returnList[[expName]][["execTime"]][["LoadingFile_From_Binary"]]=system.time(tryCatch(load(previousReturn[[expName]][["binaryObjectFileName"]]),error=function(e){cat("\n Error while reading binary file... Will read the original file instead...")returnList[[expName]][["execTime"]][["LoadingFile"]]=system.time((alignedDataList<<-readAlignedFiles(INPUTFilesList[[expName]]$folderName, INPUTFilesList[[expName]]$fileName, fileType=INPUTFilesList[[expName]]$fileType, pairedEnd=pairedEnd[[expName]], dontSort=dontSort[[expName]])))}))
#        }
#        else
#        {
#            cat("\n From file")
#            # Read the current experiment file
#            returnList[[expName]][["execTime"]][["LoadingFile"]]=system.time((alignedDataList<<-readAlignedFiles(INPUTFilesList[[expName]]$folderName, INPUTFilesList[[expName]]$fileName, fileType=INPUTFilesList[[expName]]$fileType, pairedEnd=pairedEnd[[expName]], dontSort=dontSort[[expName]])))
#        }
        
        
        returnList[[expName]][["execTime"]][[paste(programProgressName,"Reading File(s)",sep=" | ")]] <- system.time((alignedDataObject <- readAlignedData(INPUTFilesList[[expName]]$folderName, INPUTFilesList[[expName]]$fileName, fileType=INPUTFilesList[[expName]]$fileType, pairedEnds=INPUTFilesList[[expName]]$pairedEnds)))
        
        if(is.null(alignedDataObject)) stop("An error occured while reading data to process...")
        
        # Get the number of reads in the experiment
        nbReads <- length(alignedDataObject)
        if(nbReads==0) 
        {
            warningMessage="The aligned file couldn't be read properly or did not carry information about any valid aligned reads, please check your alignment statistics. Skipping..."
            cat("\nWARNING : ", warningMessage);
            warning(warningMessage);
            next
        }
        
        # Print content summary of AlignedData object (the warning is not required since the pipeline takes care of the oprhans later)
        resNULL <-  .summarizeAlignedData(alignedDataObject, noWarning=TRUE)
        
        
        # Saving the loaded data to binary R object for faster reloading in eventual recomputation
#        if("file" %in% tolower(saveBINARY))
#        {
#            binaryObjectFileName=paste(reportFilesFolder,expName, "_binaryObject.RDATA",sep="")
#            save(alignedData, file=binaryObjectFileName)
#            returnList[[expName]][["binaryObjectFileName"]]=binaryObjectFileName
#        }
#        if("return" %in% tolower(saveBINARY))
#        {
#            returnList[[expName]][["binaryObject"]]=alignedData
#        }
        
        
        ## In case a multiLocFilesList is supplied for this experiment, try to load it
        chrSplit_multiLocDataObject <- NULL
        if(expName %in% names(multiLocFilesList))
        {
            cat("\n\n Reading reads that aligned in multiple locations...")
            multiLocDataObject <- .readMultipleAlignedData(fileName=multiLocFilesList[[expName]])
            cat("\n Nb of reads aligned in multiple location :", length(multiLocDataObject))
            chrSplit_multiLocDataObject <- split(multiLocDataObject, seqnames(multiLocDataObject), drop=TRUE)
        }
        else
        {
            cat("\n No file specified for reads with multiple alignments.")
        }
        
        # Removing the chromosomes not needed
        if(nchar(removeChrNamesContaining)>0)
        {
            alignedDataObject <- dropChromosomePattern(alignedDataObject, removeChrNamesContaining, quiet=FALSE)
        }
        # Get the number of reads left after filtering
        nbReadsLeft <- length(alignedDataObject)
        if(nbReadsLeft==0) 
        {
            warningMessage="No reads left in current experiment after filtering chromosomes. Skipping..."
            cat("\nWARNING : ", warningMessage);
            warning(warningMessage);
            next
        }
        
        
        # Cleaning the chromosome name according to the prefix and suffix specified
        alignedDataObject <- normalizeChrNames(alignedDataObject, chrPrefix=INPUTFilesList[[expName]]$chrPrefix, chrSuffix=INPUTFilesList[[expName]]$chrSuffix)
        
        # Removing eventual reads with undefined strand value
        alignedDataObject <- dropUndefinedStrand(alignedDataObject, quiet=FALSE)
        nbReadsLeft <- length(alignedDataObject)
        if(nbReadsLeft==0) 
        {
            warningMessage="No reads left in current experiment after removing undefined strand values. Skipping..."
            cat("\nWARNING : ", warningMessage);
            warning(warningMessage);
            next
        }
        
        # Isolate orphans reads from paired-ends experiments
        # These reads are stored and will be eventually used in the rehab module
        chrSplit_orphanReads <- NULL
        
        # Do some specific consistency checking and organization for paired-ends experiments
        if(pairedEnds(alignedDataObject))
        {
  
            # Apply a filter that handles reads with unobvious flags/properties
            alignedDataObject <- cleanNonSimplePairs(alignedDataObject, quiet=FALSE)
            
            nbReadsLeft <- length(alignedDataObject)
            if(nbReadsLeft==0) 
            {
                warningMessage="No reads left in current experiment after cleaning pairs information. Skipping..."
                cat("\nWARNING : ", warningMessage);
                warning(warningMessage);
                next
            }
            
            
            # Isolate the orphan reads (ones with no mate)
            orphansReadsIndex <- getOrphansIndexes(alignedDataObject, quiet=FALSE)
            cat("\n     Number of orphan reads (no aligned mates) :", sum(orphansReadsIndex), "->", format((sum(orphansReadsIndex)/nbReads)*100, scientific=FALSE, digits=3), "% of aligned reads")
            if(sum(orphansReadsIndex)>0)
            {
                orphanReads <- alignedDataObject[orphansReadsIndex]
                pairedEnds(orphanReads) <- FALSE # not pairs anymore
                chrSplit_orphanReads <- split(orphanReads, seqnames(orphanReads), drop=TRUE)
                rm(orphanReads)
                gc()
                alignedDataObject <- alignedDataObject[!orphansReadsIndex]
            }
            
            nbReadsLeft <- length(alignedDataObject)
            if(nbReadsLeft==0) 
            {
                warningMessage="No reads left in current experiment after removing orphan reads. Skipping..."
                cat("\nWARNING : ", warningMessage);
                warning(warningMessage);
                next
            }
            
            
            # Filtering "outliers" inserts by size specified
            if(is.numeric(ignoreInsertsOver) && (length(ignoreInsertsOver)==1) && (ignoreInsertsOver>0))
            {
                cat("\n Filtering outliers inserts (size > ",ignoreInsertsOver,")... ", sep="")
                alignedDataObject <- filterInsertSize(alignedDataObject, 0, ignoreInsertsOver)
                cat("Done !")
            }
            
            nbReadsLeft <- length(alignedDataObject)
            if(nbReadsLeft==0) 
            {
                warningMessage="No reads left in current experiment after removing long inserts (see argument 'ignoreInsertsOver'). Skipping..."
                cat("\nWARNING : ", warningMessage);
                warning(warningMessage);
                next
            }
            
            
            # Check pairs consistency
            if(!checkPairsOK(alignedDataObject)) 
            {
                cat("\n Paired-ends reads don't seem properly sorted, trying to reorganize the data...")
                alignedDataObject <- sortByPairs(alignedDataObject, quiet=FALSE)
                if(!checkPairsOK(alignedDataObject)) stop("The dataset could not be read/sorted properly")
            }
        }
        
        cat("\n Nb of reads after filters (orphans, chromosome names, outliers) :", length(alignedDataObject), "->", format((length(alignedDataObject)/nbReads)*100, scientific=FALSE, digits=3), "% of aligned reads")
        
        #### RANGE SELECTION
        
        # Keep a reminder to know whether the ranges were estimated automatically (then we need to include the lower value of the first range) or specified manually (the user should take care of including all values in the ranges sonce lower bound is excluded)
        automaticRanges <- FALSE
        
        # In case rangeSelection is an integer, make 'rangeSelection' equal groups of size (inserts or reads if respectively pairedEnds or not)
        if(is.numeric(rangeSelection[[expName]]))
        {
            automaticRanges <- TRUE
            
            if(rangeSelection[[expName]]<2) # Should not happensince the paramter is checked to be >1 at the beginning
            {
                rangeSelection[[expName]] <- IRanges(0,-1)
            } else
            {
                cat("\n\n Preparing", rangeSelection[[expName]], "equally sized groups of", ifelse(pairedEnds(alignedDataObject),"inserts","reads"))
                groupValues <- if(pairedEnds(alignedDataObject)) isize(alignedDataObject)[isize(alignedDataObject)>0] else qwidth(alignedDataObject)
                quantRes <- quantile(groupValues, prob=seq(0,1, length.out=rangeSelection[[expName]]+1)) # compute the group limits
                # Extract the group intervals from cut function (the factor helps to remove groups eventually unused because more groups than values)
                groupsRanges <- levels(factor(cut(groupValues, breaks=unique(quantRes), include.lowest=TRUE)))
                if(length(groupsRanges)!=rangeSelection[[expName]]) cat("\n The required number of reads range based on reads/inserts size could not be met. This can happen if the user asks for too many ranges as compared to then number of different values.")
                if(length(groupsRanges)==1)
                {
                    cat("\n A single range contains all the values.")
                    rangeSelection[[expName]] <- IRanges(0,-1)
                } else
                {
                    starts <- as.numeric(sub("[\\(\\[](.+),.*", "\\1", groupsRanges))
                    ends <- as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", groupsRanges))
                    # Storing ranges
                    #rangeSelection[[expName]] <- c(IRanges(0,-1),IRanges(start=c(quantRes[1]-1,quantRes[2:(length(quantRes)-1)]), end=quantRes[2:length(quantRes)])) # Make the actual groups
                    rangeSelection[[expName]] <- c(IRanges(start=starts, end=ends), IRanges(0,-1)) # Store the groups values as IRanges used later in the code, add a "total" group
                }
            }
        }
        
        
        if(any(duplicated(rangeSelection[[expName]])))
        {
            cat("\nIgnoring one or more duplicated range of insert (or read) size. Maybe the automatic cutting has been asked for too many groups considering the values or the user defined duplicated groups in the options.")
            rangeSelection[[expName]] <- unique(rangeSelection[[expName]])
        }
        
        # Plot a summary of ranges selection on the total distribution (ignore if only total, ie. empty range)
        if((length(rangeSelection[[expName]])>1) || (!isEmpty(rangeSelection[[expName]])))
        {
            
            sizeToPlot <- if(pairedEnds(alignedDataObject)) isize(alignedDataObject)[isize(alignedDataObject)>=0] else qwidth(alignedDataObject)
            
            # Remove empty ranges (that would overlap too much with other ones)
            rangesToPlot <- rangeSelection[[expName]][!isEmpty(rangeSelection[[expName]])]
            
            colorChart <- rainbow(length(rangesToPlot))
            
            pdf(file=file.path(resultFolder, paste(expName, "_size_ranges.pdf",sep="")), width=7, height=7)
            .plotDensityRanges(sizeToPlot, ranges=as.list(as.data.frame(rbind(start(rangesToPlot), end(rangesToPlot)))), ranges.col=colorChart, labels.col="black", include.lower=c(TRUE,rep(FALSE,length(rangesToPlot)-1)), include.upper=TRUE, sep.col="darkgrey", sep.lty=2, sep.lwd=0.5, ylab="", xlab=paste("Size of ", if(pairedEnds(alignedDataObject)) "inserts" else "reads", sep=""), main="Ranges of size selection")
            dev.off()
        }
        
        # Everything is loaded and formatted, before starting processing let's allow for a selection of reads/pairs according to a size range
        # The range selection ONLY occur on reads of the experiment, MULTIREADS ARE NOT FILTERED (mainly because paired-ends multiread is not supported) !!
        for(rangeIndex in 1:length(rangeSelection[[expName]]))
        {
            
            # Iranges are not directly "loopable"...
            currentRange <- rangeSelection[[expName]][rangeIndex]
            
            # Creating a human readable range name
            rangeName <- ifelse(isEmpty(currentRange),"AllReads",paste("Range",start(currentRange),"-",end(currentRange), sep=""))
            
            cat("\n\n--------------- Range of",ifelse(INPUTFilesList[[expName]]$pairedEnds,"inserts","reads"),"size selected :", rangeName, "\n")
            
            # Saving program progress accross loops
            programProgressName_currentRange <- paste(programProgressName, rangeName, sep=" - ")
            
            # Preparing the folders for results and for figures in different ranges ("AllReads" in case of no selection)
            resultFolder_currentRange <- file.path(resultFolder,rangeName)
            reportFilesFolder <- file.path(resultFolder_currentRange, reportFilesSubFolder)
            
            # Create the figures folder, parent is created recursively
            .safeCreateFolder(reportFilesFolder)
            
            
            # Select reads in the current range
            
            rangeSelected_alignedDataObject <- NULL
            
            if(!isEmpty(currentRange))
            {
                if(pairedEnds(alignedDataObject))
                {
                    # Filter the reads (include the lower value when automatic generation of groups and treating the first one)
                    rangeSelected_alignedDataObject <- filterInsertSize(alignedDataObject, start(currentRange), end(currentRange), includeLower=(automaticRanges && (rangeIndex==1)), quiet=FALSE)
                    cat("\n Nb of reads after selection of insert range size : ", length(rangeSelected_alignedDataObject), " (",length(rangeSelected_alignedDataObject)/2, " pairs)", sep="")
                    
                    if(length(rangeSelected_alignedDataObject)>100)
                    {
                        # Plot the inserts size disctribution and the selected range
                        pdf(file=file.path(reportFilesFolder, paste(expName, "_initial_inserts_size_distribution.pdf",sep="")), width=7, height=7)
                        .plotDensityRanges(isize(alignedDataObject)[isize(alignedDataObject)>=0], ranges=c(start(currentRange), end(currentRange)), ranges.col="blue", labels.col="black", include.lower=(automaticRanges && (rangeIndex==1)), include.upper=TRUE, sep.col="darkgrey", sep.lty=2, sep.lwd=0.5, main=paste("Inserts size distribution (selection : ",length(rangeSelected_alignedDataObject)," reads)", sep=""), xlab="Inserts size (bp)", ylab="")
                        dev.off()
                    }
                    else
                    {
                        cat("\nNot enough reads to plot insert size distribution")
                    }
                }
                else # Single end, filter the READS (not the inserts) according to their size
                {
                    # Filter the reads (include the lower value when automatic generation of groups and treating the first one)
                    rangeSelected_alignedDataObject <- filterReadSize(alignedDataObject, start(currentRange), end(currentRange), includeLower=(automaticRanges && (rangeIndex==1)), quiet=FALSE)
                    cat("\n Nb of reads after selection of reads size :", length(rangeSelected_alignedDataObject))
                    
                    if(length(rangeSelected_alignedDataObject)>100)
                    {
                        # Plot the reads size disctribution
                        pdf(file=file.path(reportFilesFolder, paste(expName, "_initial_reads_size_distribution.pdf",sep="")), width=7, height=7)
                        .plotDensityRanges(qwidth(alignedDataObject), ranges=c(start(currentRange), end(currentRange)), ranges.col="blue", labels.col="black", include.lower=(automaticRanges && (rangeIndex==1)), include.upper=TRUE, sep.col="darkgrey", sep.lty=2, sep.lwd=0.5, main=paste("Reads size distribution (selection : ",length(rangeSelected_alignedDataObject)," reads)", sep=""), xlab="Reads size (bp)", ylab="")
                        dev.off()
                    }
                    else
                    {
                        cat("\nNot enough reads to plot read size distribution")
                    }
                    
                }
                
                if(length(rangeSelected_alignedDataObject)==0)
                {
                    warningMessage="No reads left after selection of current size range. Skipping..."
                    cat("\nWARNING : ", warningMessage);
                    warning(warningMessage);
                    next
                }
                
                # Plot the limits of the range selection
                #abline(v=c(start(currentRange), end(currentRange)), col="red", lty=2, lwd=2)
                
            }
            else # rangeSelected is empty, make a group with all reads
            {
                rangeSelected_alignedDataObject <- alignedDataObject
            }
            
            # Should never happen
            if(is.null(rangeSelected_alignedDataObject)) 
            {
                warningMessage="Selection of reads on size range failed. Skipping..."
                cat("\nWARNING : ", warningMessage);
                warning(warningMessage);
                next
            }
            
            
            ### If some GFF annotation files are available for this experiment, plot stats
            if(!(is.na(annotationFilesGFF[[expName]]) || is.na(annotationGenomeFiles[[expName]])))
            {
                # If the annotationGenomeFiles is not a file, maybe it's a reference to a bundled genome file, let's try to read this
                if(!file.exists(annotationGenomeFiles[[expName]]))
                {
                    precomputedReferenceFilesFolder <- "resources"
                    foundFile <- system.file(precomputedReferenceFilesFolder, paste(annotationGenomeFiles[[expName]], ".ref", sep="") ,package="Pasha")
                    
                    if(nchar(foundFile)>0) 
                    {
                        annotationGenomeFiles[[expName]] <- foundFile[1]
                    }
                    else
                    {
                        stop("The annotationGenomeFiles parameter is not a valid file or reference.")
                    }
                }
                cat("\n Plotting reads statistics on GFF files...")
                plotReadsRepartitionAnnotations(alignedData=rangeSelected_alignedDataObject,
                        gff_names_vec <-  annotationFilesGFF[[expName]], 
                        expName <- expName,
                        pdfFileName <- file.path(reportFilesFolder, paste("ReadsAnnotations_", expName, '_', rangeName, ".pdf" ,sep="")), 
                        genomeReferenceFile <- annotationGenomeFiles[[expName]])
            }
            else
            {
                cat("\n Missing GFF annotation file or genome annotation file, ignoring plot for reads occupancy statistics.")
            }
            
            
            ### Preparing data split by chromosome for further steps
            
            # Split the reads list by chromosome
            returnList[[expName]][["execTime"]][[paste(programProgressName_currentRange, "Splitting Chromosomes",sep=" | ")]]=system.time((chrSplit_alignedDataObject=split(rangeSelected_alignedDataObject, seqnames(rangeSelected_alignedDataObject), drop=TRUE)))
            
            
            # Tries to get some memory back
            rm(rangeSelected_alignedDataObject)
            if(length(rangeSelection[[expName]])<2)
            {
                rm(alignedDataObject)
            }
            gc()
            
            
            # Get the number of reads in the experiment AFTER RANGE SELECTION (used for atefact threshold definition)
            nbReads <- sum(sapply(chrSplit_alignedDataObject, length))
            
            # Get the average size of reads on ALL genome
            # Will be used for multiread pileup as the multiread file format can ommit this information. In this case we will use the average read size computed here from other reads.
            # (better to have the same for all chromosomes otherwise we might have discrepancies in the elongation size retrieved and not choose the best one, because the choice is based on the most represented one)
            averageReadSize <- trunc(mean(unlist(lapply(chrSplit_alignedDataObject, qwidth))))
            cat("\n Average reads size :",averageReadSize)
            
            
            #### REMOVING ARTEFACTS (MORE THAN threshold TAGS WITH SAME COORDINATES)
            
            # Will be used to store the indexes to remove and the elongation size computed for each artefact threshold
            #returnList[[expName]][["indexesToRemove"]] <- list()
            #returnList[[expName]][["elongationEstimation"]] <- list()
            
            for(currentIncrArtefactThrEvery in incrArtefactThrEvery[[expName]])
            {
                #browser(text=paste("Entering the loop over incrArtefactThrEvery values. current :",currentIncrArtefactThrEvery))
                
                
                chrSplit_orphansFromArtefactsReads <- NULL
                cat("\n\n------------ Artefact threshold parameter :", ifelse(is.na(currentIncrArtefactThrEvery),"Not removing artefacts",currentIncrArtefactThrEvery), "\n")
                
                # The final value used for filtering after eventual computation
                thresholdArtefacts <- NULL
                
                if(!is.na(currentIncrArtefactThrEvery))
                {
                    
                    cat("\nREMOVING ARTEFACTS ArtefactParam_", currentIncrArtefactThrEvery, sep="")
                    
                    # Define the threshold that will be used to remove the artefacts according to the number of reads
                    if(currentIncrArtefactThrEvery<0) 
                    {
                        thresholdArtefacts <- (-currentIncrArtefactThrEvery)
                    } else
                    {
                        thresholdArtefacts <- trunc(nbReads/currentIncrArtefactThrEvery)
                    }
                    
                    if(thresholdArtefacts<1) thresholdArtefacts <- 1
                    
                    # Creating a human readable artefact threshold name
                    artefactThresholdName <- paste("Artefacts threshold :", thresholdArtefacts, sep=" ")
                    # Saving program progress accross loops
                    programProgressName_artefacts <- paste(programProgressName_currentRange, artefactThresholdName, sep=" - ")
                    
                    
                    
                    ### Getting the indexes of the reads that have to be removed from the splitted object
                    
                    ## Either from a previously computed
                    # saved as the parameter used to define the threshold, for backward compatibility with versions <0.5
                    #                if(currentIncrArtefactThrEvery %in% names(previousReturn[[expName]][["indexesToRemove"]]))
                    #                {
                    #                    indexesToRemove=previousReturn[[expName]][["indexesToRemove"]][[as.character(currentIncrArtefactThrEvery)]]
                    #                }
                    #                # saved as the the threshold value, since version 0.6 because the threshold can be put manually also (as negative value)
                    #                else if(thresholdArtefacts %in% names(previousReturn[[expName]][["indexesToRemove"]]))
                    #                {
                    #                    indexesToRemove=previousReturn[[expName]][["indexesToRemove"]][[as.character(thresholdArtefacts)]]
                    #                }                
                    #                else
                    #                {
                    
                    cat("\n Threshold for artefacts removal :", thresholdArtefacts, "\n")
                    
                    # Remove the artefacts in each chromosome
                    # We don't remove it directly in case of multithreading (remove=TRUE), it would imply to transfer the objects to children processes, it's too big !
                    # Instead we only retrieve the indexes and (remove them) manipulate the big object here in the main process. It avoid copying to children (thanks
                    # to 'copy-on-write' and getting back modified result). Even if it appears more complex...
                    cat("\n Searching for artefactual piles on chromosome : ")
                    
                    returnList[[expName]][["execTime"]][[paste(programProgressName_artefacts, "Get Artefacts",sep=" | ")]] <- system.time((resultsArtefacts <- mclapply(.encapsulateListElementsInEnv(chrSplit_alignedDataObject), getArtefactsIndexes, expName, thresholdToUse=thresholdArtefacts, thresholdForStats=c(1:5,10,20,50,100), reportFilesFolder)))
                    gc()
                    
                    # Something went wrong in the fork ?
                    resultsArtefacts <- .checkErrorsInFork(resultsArtefacts)
                    
                    cat("Done !\n")
                    
                    # Merging indexes to remove from positive and negative strands of each chromosome
                    #indexesToRemove <- lapply(resultsArtefacts, function(x){return(c(x[["-"]][["indexesToRemoveFromTotal"]][[as.character(thresholdArtefacts)]],
                    #                            x[["+"]][["indexesToRemoveFromTotal"]][[as.character(thresholdArtefacts)]]))})
                    # This version does not depend on the number of strands available in the result, it merges whatever it finds in indexesToRemoveFromTotal for each chromosome
                    indexesToRemove <- lapply(resultsArtefacts, function(x) {unique(unlist(lapply(x, "[[", "indexesToRemoveFromTotal")))})
                    
                    # Removing them
                    chrSplit_noArtefact_alignedDataObject <- mapply(function(x, indexesToRemove){if(length(indexesToRemove)>0) {return(x[-indexesToRemove])} else {return(x)}}, chrSplit_alignedDataObject, indexesToRemove)
                    
                    nbReadsLeft <- length(chrSplit_noArtefact_alignedDataObject)
                    if(nbReadsLeft==0) 
                    {
                        warningMessage="No reads left in current range after removing artefact reads. Skipping..."
                        cat("\nWARNING : ", warningMessage);
                        warning(warningMessage);
                        next
                    }
                    
                    
                    # Handle orphans from artefacts removal
                    if(INPUTFilesList[[expName]]$pairedEnds)
                    {
                        # Get orphans indexes created by removing artefacts
                        chrSplit_orphansIndexes <- lapply(chrSplit_noArtefact_alignedDataObject,getOrphansIndexes)
                        
                        # Check that there has been some orphans found before trying to get them, we don't want to alter the NULL value in case there is no orphan (it will be used to decide whether to compute or not the piling)
                        if(any(sapply(chrSplit_orphansIndexes, any)))
                        {
                            # Isolate the orphans for future rehabilitation
                            chrSplit_orphansFromArtefactsReads <- mapply("[", chrSplit_noArtefact_alignedDataObject, chrSplit_orphansIndexes) # The "[" operator of the class automatically return an empty object in case there is no indexes 
                            chrSplit_orphansFromArtefactsReads <- lapply(chrSplit_orphansFromArtefactsReads, "pairedEnds<-", FALSE) # Not pairs anymore
                            
                            # Removing them (orphans) from general object, IMPORTANT : chrSplit_orphansIndexes IS A LOGICA VECTOR, not indexes
                            chrSplit_noArtefact_alignedDataObject <- mapply(function(x, indexesToRemove){if(length(indexesToRemove)>0) {return(x[-indexesToRemove])} else {return(x)}}, chrSplit_noArtefact_alignedDataObject, lapply(chrSplit_orphansIndexes,which))    
                        }
                    }
                    
                    nbReadsLeft <- length(chrSplit_noArtefact_alignedDataObject)
                    if(nbReadsLeft==0) 
                    {
                        warningMessage="No reads left in current range after removing orphans generated by artefact removal. Skipping..."
                        cat("\nWARNING : ", warningMessage);
                        warning(warningMessage);
                        next
                    }

                    
                    ### TODO : plot artefacts stats by chromosome
                    
                    
                    # Saving the indexes to remove for eventual other computation on the same dataset (will help to avoid this long step)
                    #returnList[[expName]][["indexesToRemove"]][[as.character(currentIncrArtefactThrEvery)]]=indexesToRemove # version <0.6
                    #returnList[[expName]][["indexesToRemove"]][[as.character(thresholdArtefacts)]]=indexesToRemove # version >=0.6
                    
                    # Removing the indexes which fall in artefacts (classic mapply, not multithreaded, see precedent comment)
                    # Mapply will put in correspondance the elements of the first list to the ones of the second, we put indexes as negative to remove them with "["
                    #returnList[[expName]][["execTime"]][[paste("ArtefactsRemove", "_ArtefactParam_",currentIncrArtefactThrEvery, sep="")]]=system.time((chrSplittedDataMinusArtefacts=mapply(function(data, chrNames){if(length(indexesToRemove[[chrNames]])>0) return(data[-indexesToRemove[[chrNames]]]) else return(data)}, chrSplittedData, names(chrSplittedData))))
                    #returnList[[expName]][["execTime"]][[paste("ArtefactsRemove", "_ArtefactParam_",currentIncrArtefactThrEvery, sep="")]]=system.time((chrSplittedDataMinusArtefacts=mapply("[",chrSplittedData,lapply(indexesToRemove,"-"))))
                    
                    nbRemoved <- length(unlist(indexesToRemove))
                    nbRemaining <- sum(sapply(chrSplit_noArtefact_alignedDataObject,length))
                    
                    cat("\n Nb of reads removed :", nbRemoved, "->", format((nbRemoved/nbReads)*100, scientific=FALSE, digits=3), "% of aligned reads")
                    cat("\n Nb of reads remaining :", nbRemaining, "->", format((nbRemaining/nbReads)*100, scientific=FALSE, digits=3), "% of aligned reads")
                    
                } # /REMOVING ARTEFACTS
                else
                {
                    
                    # Saving program progress accross loops
                    programProgressName_artefacts <- paste(programProgressName_currentRange, "No Artefact Removal", sep=" - ")
                    
                    # Save unchanged object
                    chrSplit_noArtefact_alignedDataObject <- chrSplit_alignedDataObject
                }
                
                if(length(incrArtefactThrEvery[[expName]])<2)
                {
                    rm(chrSplit_alignedDataObject)
                }
                gc()
                
                #### ELONGATION SIZE ESTIMATION
                
                for(currentElongationSize in elongationSize[[expName]])
                {
                    cat("\n\n--------- Elongation/Shifting :", ifelse(is.na(currentElongationSize),ifelse(INPUTFilesList[[expName]]$pairedEnds,"based on pairs","automatic estimation"),currentElongationSize),"\n")
                    
                    elongationName <- NULL
                    
                    # Flag to mark whether the elongation comes from a manual parameter or by estimation (in case estimation is the same as a manual specification), used to generate output files name
                    manualElongation <- TRUE
                                        
                    if(is.na(currentElongationSize))
                    {
                        manualElongation <- FALSE
                        if(!INPUTFilesList[[expName]]$pairedEnds)
                        {
                            cat("\nELONGATION SIZE ESTIMATION")
                            
                            # In case the size has been estimated in a previous computation, retrieve it from the 'previousReturn' object
                            #                    if(currentIncrArtefactThrEvery %in% names(previousReturn[[expName]][["elongationEstimation"]]))
                            #                    {
                            #                        currentElongationSize <- previousReturn[[expName]][["elongationEstimation"]][[as.character(currentIncrArtefactThrEvery)]]
                            #                    }
                            #                    else
                            #                    {
                            
                            cat("\n Estimating elongation/shifting size on chromosome : ")
                            
                            # Estimate the elongation size for each chromosome
                            returnList[[expName]][["execTime"]][[paste(programProgressName_artefacts, "Elongation/Shifting Estimation",sep=" | ")]]=system.time((elongationEstimationList=mclapply(.encapsulateListElementsInEnv(chrSplit_noArtefact_alignedDataObject), estimateElongationSize, expName, elongationEstimationRange["mini"], elongationEstimationRange["maxi"], elongationEstimationRange["by"], averageReadSize, reportFilesFolder)))
                            gc()
                            
                            # Something went wrong in the fork ?
                            elongationEstimationList <- .checkErrorsInFork(elongationEstimationList)
                            
                            cat("Done !\n")
                            
                            
                            # Get the result as a named vector (avoid the redundant name coming from the named vector returned by the function and from the named list from apply)
                            # Get the names from the list object so all chromosomes with no reads won't be considered anymore (NULL results are automatically removed by unlist)
                            elongationEstimations <- unlist(lapply(elongationEstimationList, unname))
                            
                            # Get the number of each score
                            elongationCounts <- table(elongationEstimations)
                            # Identify the elongations that are less than twice the read size (they are considered as suspicious values and will be weighted down for final decision)
                            smallElongationIndex <- (as.numeric(names(elongationCounts)) < (averageReadSize*2))
                            # Weight down the counts for small elongation values (they have to be represented twice as much as other scores for being selected)
                            elongationCounts[smallElongationIndex] <- elongationCounts[smallElongationIndex]/2;
                            
                            # Get the most represented size among chromosomes
                            currentElongationSize <- as.numeric(names(which.max(elongationCounts)))
                            
                            if(currentElongationSize < (averageReadSize*2))
                            {
                                warningMessage <- "Small elongation values (less than twice the reads size) were the most represented among all chromosomes (even after these values were weighted down for selection). This can reflect anomalous fragment size estimation, please consider checking elongation reports for individual chromosomes and eventually specify a manual value..."
                                # Write warning in imediate log file
                                cat("\n WARNING :", warningMessage)
                                # Report warning programatically
                                warning(warningMessage)
                            }
                            
                            cat("\n Estimation summary (bp) :\n    ")
                            
                            summaryText <- c(paste(names(elongationEstimations), " --> ", elongationEstimations, sep=""), paste("----> Most represented :", currentElongationSize,"bp"))
                            
                            cat(summaryText, sep="\n    ")
                            
                            # Write the summary text file for all chromosomes
                            writeLines(con=file.path(reportFilesFolder,"Shifting_Estimation_Summary.txt"),text=summaryText)
                            #}
                            
                            # Saving the result for eventual other computation on the same dataset (will help to avoid this long step)
                            #returnList[[expName]][["elongationEstimation"]][[as.character(currentIncrArtefactThrEvery)]]=currentElongationSize
                            
                            # Creating a human readable elongation/shifting size
                            elongationName <- paste("Estimated elongation/shifting", currentElongationSize, sep=" ")
                            
                            
                        } 
                        else # Paired Ends
                        {
                            # In this case the elongation size is informative (mainly here to plot a graph distribution)
                            # A decision (a number) will be taken however based on the distribution that could be useful for eventual future 'orphan reads' rehabilitation
                            currentElongationSize <- trunc(median(unlist(lapply(lapply(chrSplit_noArtefact_alignedDataObject,isize), abs))))
                            
                            cat("\n Median elongation size computed using pairs (useful for eventual elongation of orphans or multireads) :", currentElongationSize)
                            
                            
                            # Plotting
                            # Must select for dataset (chromosomes) that have enough valid data to plot a density (at least 2 values to plot the density in positive insert sizes) 
                            chrSplit_noArtefact_alignedDataObject_validValues <- mapply("[", chrSplit_noArtefact_alignedDataObject, lapply(chrSplit_noArtefact_alignedDataObject, function(x) {return(isize(x)>0)}))
                            chrSplit_noArtefact_alignedDataObject_min2values <- chrSplit_noArtefact_alignedDataObject_validValues[sapply(chrSplit_noArtefact_alignedDataObject_validValues, length)>2]
                            if(length(chrSplit_noArtefact_alignedDataObject_min2values)>0)
                            {
                                densityInsertSize <- lapply(lapply(chrSplit_noArtefact_alignedDataObject_min2values, function(x){return(isize(x)[isize(x)>0])}), density)
                                
                                ### Isolating the median of max density points
                                # Get the max density points
                                densityMaxPoints <- sapply(densityInsertSize, function(densityData){return(densityData$x[which.max(densityData$y)])})
                                
                                # Preparing merged plotting
                                max.y <- max(sapply(densityInsertSize,"[[", "y"))
                                max.x <- max(sapply(densityInsertSize,"[[", "x"))
                                
                                pdf(file=file.path(reportFilesFolder, paste(expName, "_insert_size_distribution_by_chromosome",currentElongationSize,".pdf",sep="")), width=7, height=7)    
                                # Set the window
                                plot(NULL, xlim=c(0,max.x), ylim=c(0,max.y), main=paste(expName, " insertSize distributions (",currentElongationSize,")", sep=""), xlab="Insert size", ylab="")
                                
                                # Plot the curves
                                resNULL <- lapply(densityInsertSize, function(densityData){points(densityData, col= "#55555555", lty=2, type='l')})
                                abline(v=currentElongationSize, col="red", lwd=2)
                                
                                dev.off()
                            }
                            else
                            {
                                cat("\nNot enough reads to plot insert size distribution")
                            }
                            
                            # Creating a human readable elongation/shifting size
                            elongationName <- paste("Estimated elongation/shifting (based on pairs for orphans and multireads)", currentElongationSize, sep=" ")                    
                            
                        } # /ELONGATION SIZE ESTIMATION
                        
                    }
                    else # elongation/shifting size manually specified
                    {
                        if(INPUTFilesList[[expName]]$pairedEnds && (currentElongationSize != 0))
                        {
                            cat("\n WARNING : a manual elongation size (other than 0) has been specified for a paired-ends dataset. Note that this value will only be used for eventual orphans and multireads (standard piling will use length information carried by pairs)...");
                        }
                        # Creating a human readable elongation/shifting size
                        elongationName <- paste("Manual elongation/shifting", currentElongationSize, sep=" ")
                    }
                    
                    # Saving program progress accross loops
                    programProgressName_elongation <- paste(programProgressName_artefacts, elongationName, sep=" - ")
                    
                    #### RESULT GENERATION
                    
                    # This variable will store intermediary piled-up results for each category ("mergedReads", "unireads", "multiReads", "orphans", "orphansFromArtefacts") 
                    piledRleData <- list()
                    # mergedReads
                    # unireads
                    # multiReads
                    # orphans
                    # orphansFromArtefacts
                    
                    
                    if(any(c(WIGvs,WIGfs,GFF, BED, BIGWIG)))
                    {
                        cat("\n\n COMPUTING PILED FILES\n")
                        
                        
                        
                        cat("\n Piling UNIREADS reads on chromosome : ")
                        
                        # Computing the pileup vector coming from the experiment and putting it in the uniread slot of the list
                        returnList[[expName]][["execTime"]][[paste(programProgressName_elongation, "Piling Unireads",sep=" | ")]]=system.time((piledRleData[["unireads"]]=mclapply(.encapsulateListElementsInEnv(chrSplit_noArtefact_alignedDataObject), generatePiled, ifelse(currentElongationSize==0, 0,if(INPUTFilesList[[expName]]$pairedEnds) NA else currentElongationSize), averageReadSize, INPUTFilesList[[expName]]$midPoint)))
                        if(length(elongationSize[[expName]])<2)
                        {
                            rm(chrSplit_noArtefact_alignedDataObject)
                        }
                        gc()
                        
                        # Something went wrong in the fork ?
                        piledRleData[["unireads"]] <- .checkErrorsInFork(piledRleData[["unireads"]])
                        
                        
                        cat("Done !\n")
                        
                        
                        
                        ### Preparing which piling will be done
                        
                        # Vector checking (bool) that the data is available for each eventual piling step that will occur
                        dataAvailablePiling <- c(orphans=((!is.null(chrSplit_orphanReads)) && (length(chrSplit_orphanReads)>0)), 
                                orphansFromArtefacts=((!is.null(chrSplit_orphansFromArtefactsReads)) && (length(chrSplit_orphansFromArtefactsReads)>0)), 
                                multiReads=((!is.null(chrSplit_multiLocDataObject)) && (length(chrSplit_multiLocDataObject)>0)))
                        
                        userAskingPiling <- c(orphans=("orphans" %in% rehabilitationStep), 
                                orphansFromArtefacts=("orphansFromArtefacts" %in% rehabilitationStep), 
                                multiReads=TRUE) # as soon as the user gives a file (checked by other test), generate it
                        
                        # Vector that summarize if the piling will be done or not for eache eventual step
                        doPiling <- dataAvailablePiling & userAskingPiling
                        
                        
                        ### Doing the alternative/optional pilings
                        
                        # Rehab orphans reads coming from broken pairs
                        if(doPiling["orphans"])
                        {
                            cat("\n Rehabilitation and piling of orphan reads on chromosome : ")
                            returnList[[expName]][["execTime"]][[paste(programProgressName_elongation, "Piling Orphans",sep=" | ")]]=system.time((piledRleData[["orphans"]]=mclapply(.encapsulateListElementsInEnv(chrSplit_orphanReads), generatePiled, currentElongationSize, averageReadSize, INPUTFilesList[[expName]]$midPoint)))
                            if(length(elongationSize[[expName]])<2)
                            {
                                rm(chrSplit_orphanReads)
                            }
                            gc()
                            # Something went wrong in the fork ?
                            piledRleData[["orphans"]] <- .checkErrorsInFork(piledRleData[["orphans"]])
                            cat("Done !\n")
                        }
                        
                        # Rehab orphans reads coming from broken pairs during artefact removing step
                        if(doPiling["orphansFromArtefacts"])
                        {
                            cat("\n Rehabilitation and piling of orphan reads from artefact removing step on chromosome : ")
                            returnList[[expName]][["execTime"]][[paste(programProgressName_elongation, "Piling Orphans From Artefacts",sep=" | ")]] <- system.time((piledRleData[["orphansFromArtefacts"]]=mclapply(.encapsulateListElementsInEnv(chrSplit_orphansFromArtefactsReads), generatePiled, currentElongationSize, averageReadSize, INPUTFilesList[[expName]]$midPoint)))
                            if(length(elongationSize[[expName]])<2)
                            {
                                rm(chrSplit_orphansFromArtefactsReads)
                            }
                            gc()
                            # Something went wrong in the fork ?
                            piledRleData[["orphansFromArtefacts"]] <- .checkErrorsInFork(piledRleData[["orphansFromArtefacts"]])
                            cat("Done !\n")
                        }
                        
                        # Process reads related to current experiment that aligned in multiple locations
                        if(doPiling["multiReads"])
                        {
                            cat("\n Piling multireads (reads aligned in multiple locations) on chromosome : ")
                            returnList[[expName]][["execTime"]][[paste(programProgressName_elongation, "Piling Multireads",sep=" | ")]] <- system.time((piledRleData[["multiReads"]] <- mclapply(.encapsulateListElementsInEnv(chrSplit_multiLocDataObject), generatePiled, currentElongationSize, averageReadSize, INPUTFilesList[[expName]]$midPoint)))
                            if(length(elongationSize[[expName]])<2)
                            {
                                rm(chrSplit_multiLocDataObject)
                            }
                            gc()
                            # Something went wrong in the fork ?
                            piledRleData[["multiReads"]] <- .checkErrorsInFork(piledRleData[["multiReads"]])
                            cat("Done !\n")
                        }
                        
                        # In case there has been several ones, compute the merge of all for the final file
                        if(any(doPiling))
                        {
                            cat("\n Preparing and merging final piled vector... ")
                            
                            
                            #if(length(piledRleData)>1)
                            #{
                            
                            # Get a list of all chromosomes in all piled lists (can't use mapply in case of heterogeneous chromosomes lists)
                            chrList <- unique(unlist(lapply(piledRleData, names)))
                            
                            # Trick to keep the names after lapply on the chromosomes names
                            names(chrList) <- chrList
                            #browser()
                            returnList[[expName]][["execTime"]][[paste(programProgressName_elongation, "Merging Pileups",sep=" | ")]] <- system.time((piledRleData[["mergedReads"]] <- mclapply(chrList, function(currentChr)
                                                        {
                                                            piledSet <- c(piledRleData["unireads"],piledRleData[names(doPiling)[doPiling]]) # select only the ones that were supposedly computed
                                                            piledSetChrSelection <- lapply(piledSet, "[[", currentChr) # Select the chromosome of interest
                                                            
                                                            # Removing the ones that had no information for the current chromosome (in case of very unbalanced datasets)
                                                            piledSetChrSelection <- piledSetChrSelection[!sapply(piledSetChrSelection,is.null)]
                                                            
                                                            # If there is at least one dataset covering this chromosome (should always be the case because we only treat chromosomes based on what was found in the piled object)
                                                            if(length(piledSetChrSelection)>0)
                                                            {
                                                                
                                                                # Resize to maxSize by appending 0s to all Rle vectors to avoid recycling)
                                                                maxSize <- max(sapply(piledSetChrSelection, length))
                                                                piledSetChrSelection <- lapply(piledSetChrSelection, function(x){return(append(x, rep(0,maxSize-length(x))))})
                                                                
                                                                # Create the resulting sum vector
                                                                mergedRle <- Reduce('+', piledSetChrSelection)
                                                                
                                                                # Adding the chromosome name to metadata after all computations, altering a Rle removes its metadata...
                                                                metadata(mergedRle) <- list("chr"=currentChr) # Keep the chromosome name information for further processing
                                                                
                                                                return(mergedRle)
                                                            }
                                                            else
                                                            {
                                                                # This will be automatically removed from returned results in .checkErrorsInFork
                                                                return(NULL)
                                                            }
                                                        }))) 
                            
                            # Something went wrong in the fork ?
                            piledRleData[["mergedReads"]] <- .checkErrorsInFork(piledRleData[["mergedReads"]])
                            
                            
                            if(!keepTemp)
                            {
                                piledRleData <- piledRleData["mergedReads"] # returns a sublist containing only merged dataset
                            }
                            #}
                            #else
                            #{
                            #    names(piledRleData)="mergedReads"
                            #}
                            
                            cat("Done !")
                            
                        }
                        
                        
                        
                        #### WRITING FILES 
                        
                        # Create or not a subfolder in resultFolder depending on the number of pileup categories
                        pileupSubcategoriesFolder <- file.path(resultFolder_currentRange, "PileupSubcategories")
                        if(length(piledRleData)>1)
                        {
                            .safeCreateFolder(pileupSubcategoriesFolder)
                        }
                        
                        
                        # Will loop across all piled categories to make the final files
                        for(currentPileupCategory in names(piledRleData))
                        {
                            cat("\n\n------ Writing on disk :", currentPileupCategory,"\n")
                            
                            # Saving program progress accross loops
                            programProgressName_pileup <- paste(programProgressName_elongation, currentPileupCategory, sep=" - ")
                            
                            
                            
                            # Put in the subfolder if we deal with subcategories
                            finalResultFolder <- ifelse((currentPileupCategory=="mergedReads") || (length(piledRleData)==1), resultFolder_currentRange, pileupSubcategoriesFolder)
                            
                            
                            # Create a base filename that will be used to create the final filenames and the track names in WIGs
                            #baseFileName <- paste(expName, ifelse(length(piledRleData)>1, paste("_", currentPileupCategory,sep=""), ""), ifelse(INPUTFilesList[[expName]]$midPoint, "_sh", "_el"), currentElongationSize,ifelse(!is.na(currentIncrArtefactThrEvery),paste("_AThr",thresholdArtefacts, sep=""),""),ifelse(INPUTFilesList[[expName]]$midPoint, "_MIDPOINT", ""), sep="")
                            baseFileName <- paste(expName, "_", currentPileupCategory, ifelse(INPUTFilesList[[expName]]$midPoint, "_sh", "_el"), ifelse(INPUTFilesList[[expName]]$pairedEnds, "PairsAnd", ""),ifelse(manualElongation, "Manual", "Est"),currentElongationSize,ifelse(!is.na(currentIncrArtefactThrEvery),paste("_AThr",thresholdArtefacts, sep=""),""),ifelse(INPUTFilesList[[expName]]$midPoint, "_MIDPOINT", ""), sep="")
                            
                            # Generation of the WIG variable step which IS NOT dependant of the binSize
                            if(WIGvs || BIGWIG)
                            {
                                cat("\n Writing WIG variable steps for chromosome : ")

                                returnList[[expName]][["execTime"]][[paste(programProgressName_pileup, "Writing WIGvs chromosomes",sep=" | ")]]=system.time((resultingFiles_WIGvs=mclapply(piledRleData[[currentPileupCategory]], .writeWIGvs_chr, baseFileName, finalResultFolder, compatibilityOutputWIG)))
                                gc()
                                
                                # Something went wrong in the fork ?
                                resultingFiles_WIGvs <- .checkErrorsInFork(resultingFiles_WIGvs)
                                
                                cat("Done !\n")
                                
                                
                                cat("\n Merging WIG variable step temporary result files")
                                
                                # Get the names of wig files for chromosomes
                                chromosomesFiles <- mixedsort(unname(unlist(resultingFiles_WIGvs)))
                                
                                # Create the merged result File
                                resultFileName <- file.path(finalResultFolder, paste("WIGvs_", baseFileName, ".wig", sep=""))
                                
                                # Clean an eventual previous file with the same name
                                if(file.exists(resultFileName))
                                {
                                    file.remove(resultFileName)
                                }
                                
                                # Write the track line in first line of the result file (chromosomes files will be appended after it, when compatibility mode on this is written by writeWIG function for each chromosome) 
                                if(!compatibilityOutputWIG)
                                {
                                    fileCon <- file(resultFileName, open="wb")
                                    writeLines(paste("track type=wiggle_0 name=\"", baseFileName, "\"", sep=""), con=fileCon, sep="\n")
                                    close(fileCon)
                                }
                                
                                # Concatenate the variable step wig file (thanks to rle encoding)
                                returnList[[expName]][["execTime"]][[paste(programProgressName_pileup, "Appending WIGvs Individual Chromosomes",sep=" | ")]] <- system.time(file.append(resultFileName, chromosomesFiles))
                                
                                #                        if(!keepTemp)
                                #                        {
                                # Erase the temporary (chromosomes splitted) WIG files
                                file.remove(chromosomesFiles)
                                #                        }
	
								if(BIGWIG)
								{ 
                                    
                                    
                                    # Compute the minimal chromosome length from piled scores in order to generate the seqinfo object required for bigwig conversion
                                    chromLengths=sapply(lapply(piledRleData[[currentPileupCategory]], runLength), sum);
                                    seqinfo_object=Seqinfo(seqnames=names(chromLengths), seqlengths=chromLengths, isCircular=rep(FALSE, length(chromLengths)), genome="genome")
                                    
									cat("\n Writing BIGWIG file : ")	
									wigToBigWig(resultFileName, seqinfo_object, gsub("\\.wig$", ".bw", resultFileName));
									
									if(!WIGvs)
									{
										file.remove(resultFileName);
									}
								}
                            }
                            
                            # Generation of the BEDGRAPH which IS NOT dependant of the binSize
                            if(BED)
                            {
                                cat("\n Writing Bedgraph for chromosome : ")
                                
                                returnList[[expName]][["execTime"]][[paste(programProgressName_pileup, "Writing Bedgraph chromosomes",sep=" | ")]] <- system.time((resultingFiles_BED <- mclapply(piledRleData[[currentPileupCategory]], .writeBED_chr, baseFileName, finalResultFolder)))
                                gc()
                                
                                # Something went wrong in the fork ?
                                resultingFiles_BED <- .checkErrorsInFork(resultingFiles_BED)
                                
                                cat("Done !\n")
                                
                                
                                cat("\n Merging Bedgraph temporary result files")
                                
                                # Get the names of wig files for chromosomes
                                chromosomesFiles <- mixedsort(unname(unlist(resultingFiles_BED)))
                                
                                # Create the merged result File
                                resultFileName <- file.path(finalResultFolder, paste("BED_", baseFileName, ".bed", sep=""))
                                
                                # Clean an eventual previous file with the same name
                                if(file.exists(resultFileName))
                                {
                                    file.remove(resultFileName)
                                }
                                
                                # Write the track line in first line of the result file
                                fileCon <- file(resultFileName, open="wb")
                                writeLines(paste("track type=bedGraph name=\"", baseFileName, "\"", sep=""), con=fileCon, sep="\n")
                                close(fileCon)
                                
                                # Concatenate the GFF file (for peak detection with CoCas)
                                returnList[[expName]][["execTime"]][[paste(programProgressName_pileup, "Appending Bedgraph Individual Chromosomes",sep=" | ")]] <- system.time(file.append(resultFileName, chromosomesFiles))
                                
                                #                                if(!keepTemp)
                                #                                {
                                # Erase the temporary (chromosomes splitted) GFF files
                                file.remove(chromosomesFiles)
                                #                                }
                            }
                            
                            
                            # Loop on the eventually different binsizes specified by user
                            for(currentBinSize in binSize[[expName]])
                            {
                                cat("\n\n--- Bin Size :", currentBinSize,"\n")
                                
                                # Creating a human readable Binsize
                                binSizeName <- paste("Binsize", currentBinSize, sep=" ")
                                # Saving program progress accross loops
                                programProgressName_binSize <- paste(programProgressName_pileup, binSizeName, sep=" - ")
                                
                                
                                # Generation of the WIG FS and the GFF which ARE dependant of the binSize
                                if(any(c(WIGfs,GFF)) & (!is.na(currentBinSize)))
                                {
                                    
                                    baseFileNameBinned <- paste(baseFileName,"_bin",currentBinSize,sep="")
                                    
                                    ## BINNING
                                    
                                    
                                    binningVec <- function(piledRleData, binSize)
                                    {
                                        currentChr <- metadata(piledRleData)[["chr"]]
                                        
                                        cat(paste(currentChr, "... ", sep=""))
                                        res <- binVector(as.numeric(piledRleData), binSize)
                                        res <- Rle(res)
                                        metadata(res) <- list("chr"=currentChr)
                                        return(res)
                                    }
									
									
									
                                    # We keep Rle data representation for memory saving despite the fact that we will have to convert it for writing
                                    # Moreover we'll be able to use metadata to store the chromosome name
                                    if(currentBinSize>1)
                                    {
                                        cat("\n Binning chromosome (",currentBinSize,"bp) : ",sep="")
                                        
                                        returnList[[expName]][["execTime"]][[paste(programProgressName_binSize, "Binning",sep=" | ")]] <- system.time((binnedDataRle <- mclapply(piledRleData[[currentPileupCategory]], binningVec,  currentBinSize)))
                                        gc()
                                        
                                        # Something went wrong in the fork ?
                                        piledRleData[[currentPileupCategory]] <- .checkErrorsInFork(piledRleData[[currentPileupCategory]])
                                        
                                        cat("Done !\n")
                                    }
                                    else
                                    {
                                        binnedDataRle <- piledRleData[[currentPileupCategory]]
                                    }
                                    
                                    
                                    
                                    if(WIGfs)
                                    {
                                        
                                        cat("\n Writing WIG fixed steps (bins) for chromosome : ")
                                        
                                        returnList[[expName]][["execTime"]][[paste(programProgressName_binSize, "Writing WIGfs chromosomes",sep=" | ")]] <- system.time((resultingFiles_WIGfs <- mclapply(binnedDataRle, .writeWIGfs_chr, baseFileNameBinned, currentBinSize, finalResultFolder, compatibilityOutputWIG)))
                                        gc()
                                        
                                        # Something went wrong in the fork ?
                                        resultingFiles_WIGfs <- .checkErrorsInFork(resultingFiles_WIGfs)

                                        cat("Done !\n")
                                        
                                        
                                        cat("\n Merging WIG fixed step temporary result files")
                                        
                                        # Get the names of wig files for individual chromosomes (and sort them by chromosome name)
                                        chromosomesFiles <- mixedsort(unname(unlist(resultingFiles_WIGfs)))
                                        
                                        # Create the merged result File
                                        resultFileName <- file.path(finalResultFolder, paste("WIGfs_", baseFileNameBinned, ".wig", sep=""))
                                        
                                        # Clean an eventual previous file with the same name
                                        if(file.exists(resultFileName))
                                        {
                                            file.remove(resultFileName)
                                        }
                                        
                                        # Write the track line in first line of the result file (chromosomes files will be appended after it, when compatibility mode on this is written by writeWIG function for each chromosome) 
                                        if(!compatibilityOutputWIG)
                                        {
                                            fileCon <- file(resultFileName, open="wb")
                                            writeLines(paste("track type=wiggle_0 name=\"", baseFileNameBinned, "\"", sep=""), con=fileCon, sep="\n")
                                            close(fileCon)
                                        }
                                        
                                        
                                        # Concatenate the fixed step (currentBinSize) wig files
                                        returnList[[expName]][["execTime"]][[paste(programProgressName_binSize, "Appending WIGfs Individual Chromosomes",sep=" | ")]] <- system.time(file.append(resultFileName, chromosomesFiles))
                                        
                                        #                                if(!keepTemp)
                                        #                                {
                                        # Erase the temporary (chromosomes splitted) WIG files
                                        file.remove(chromosomesFiles)
                                        #                                }
                                    }
                                    
                                    
                                    if(GFF)
                                    {
                                        cat("\n Writing GFF for chromosome : ")
                                        
                                        returnList[[expName]][["execTime"]][[paste(programProgressName_binSize, "Writing GFF chromosomes",sep=" | ")]] <- system.time((resultingFiles_GFF <- mclapply(binnedDataRle, .writeGFF_chr, baseFileNameBinned, currentBinSize, finalResultFolder)))
                                        gc()
                                        
                                        # Something went wrong in the fork ?
                                        resultingFiles_GFF <- .checkErrorsInFork(resultingFiles_GFF)
                                        
                                        cat("Done !\n")
                                        
                                        
                                        cat("\n Merging GFF temporary result files")
                                        
                                        # Get the names of wig files for chromosomes
                                        chromosomesFiles <- mixedsort(unname(unlist(resultingFiles_GFF)))
                                        
                                        # Create the merged result File
                                        resultFileName <- file.path(finalResultFolder, paste("GFF_", baseFileNameBinned, ".gff", sep=""))
                                        
                                        # Clean an eventual previous file with the same name
                                        if(file.exists(resultFileName))
                                        {
                                            file.remove(resultFileName)
                                        }
                                        
                                        # Concatenate the GFF file (for peak detection with CoCas)
                                        returnList[[expName]][["execTime"]][[paste(programProgressName_binSize, "Appending GFF Individual Chromosomes",sep=" | ")]] <- system.time(file.append(resultFileName, chromosomesFiles))
                                        
                                        #                                if(!keepTemp)
                                        #                                {
                                        # Erase the temporary (chromosomes splitted) GFF files
                                        file.remove(chromosomesFiles)
                                        #                                }
                                    }

                                    
                                } # /WIG FS AND GFF GENERATION
                                
                            } # /binSize
                        } # currentPileupCategory
                        
                    } # any(c(WIGvs,WIGfs,GFF))
                    
                } # /ElongationSize
                
            } # /ArtefactsThreshold
            
        } # /currentRange
        
        ts <- paste("##------", date(), "------##")
        cat("\n",ts,"\n")
        
        .printTimes(returnList[expName])
        
        # Releasing the log file for this experiment (sink automatically manages it like a fifo)
        sink()
        
    } # /expName
    
    
    
    .printTimes(returnList)
    
    ts <- paste("##------", date(), "------##")
    cat("\n\n",ts,"\n")
    
    cat("\n\nDone, successfully processed",length(INPUTFilesList),"experiment(s).")
    
    
    # Releasing the log files
    if(generalLogFile)
    {
        sink()
    }
    
    return(returnList)
    
} # /processBAM


