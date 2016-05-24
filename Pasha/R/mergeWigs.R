###############################
# ChipSeq replicate merging
# Merging wig files
# 20/11/2008 - R.Fenouil
###############################


mergeWigs = function(files, binSize=50, outputFolder="./") ## mean, sum   ## subtract, ratio, remove
{
    
#Arguments validation :
# - ip non empty
# - files exists
# - files is named
    
    dir.create(outputFolder, recursive=TRUE)
    
    # Computing separate biological replicates (or experiments) => one output file per element
    for(expName in names(files))
    {
        title <- paste("------ EXP ",expName," ------",sep="")
        cat("\n")
        cat(paste(rep("-",nchar(title)),collapse=""),"\n")
        cat(title,"\n")
        cat(paste(rep("-",nchar(title)),collapse=""),"\n")
        
        cat("\n  Reading IPs :\n")
        
        ips <- list()
        
        #reading files and get it as a list of list (one experiment per element then one chromosome per element)
        for(index in 1:length(files[[expName]]))
        {
            cat("    Reading ", files[[expName]][index], " ... ")
            ips[[index]] <- readWIG(files[[expName]][index])
            cat("Done.\n")
        }
        
        chrNamesList <- lapply(ips, names)
        
        # Get the chromosomes names in common for all experiments to merge
#        chrNames=do.call(intersect, chrNamesList)
        chrNames <- unique(unname(unlist(chrNamesList)))
        for(i in 1:length(chrNamesList))
        {
            chrNames <- intersect(chrNames, chrNamesList[[i]])
        }
        
        # Ordering the chromosomes for all the experiments
        ips <- lapply(ips, "[", chrNames)
        
        
        
        ipRes <- list()
        
        cat("\n  Merging IPs :\n")
        
        
        # computing merge for each chromosome
        for(chrCurrent in chrNames)
        {
            
            cat("    Merging chromosome ",chrCurrent," ... ")
            
            # Isolate the current chromosome among all experiments
            chrDataList <- lapply(ips, "[[", chrCurrent)
            
            # Get the max size for this chromosome among the experiments (can be different because of the pileup)
            maxCurrentSize <- max(sapply(chrDataList, length))
            
            # Setting the size of all the chromosomes to the longer one (padding with 0s)
            chrDataList <- lapply(chrDataList, function(x,y)
                    {
                        length(x) <- y
                        x[is.na(x)] <- 0
                        return(x)
                    }, maxCurrentSize)
            
            # Putting all the values as columns of a matrix
            chrDataList <- do.call(cbind,chrDataList)
            
            # applying the merging method to the matrix lines in order to retrieve a vector 
            # corresponding to the result of merging. Putting this vector in the result wigList.
#            ipRes[[chrCurrent]]=apply(chrDataList,1,mergingMethod)
            ipRes[[chrCurrent]] <- rowMeans(chrDataList)
            
            cat("Done.\n")
        }
        
        # we are finished with the IPs, let's go for inputs
        rm(ips)
        invisible(gc())
        
        
        cat("\n  Writing WIG result file ... ")
        
        writeWIG(ipRes, expName, folder=outputFolder, fixedStep=binSize)
        
        cat("Done.\n")
        
    } # for experiments
    
    return(NULL)
}



