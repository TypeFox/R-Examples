
multiread_RemoveArtifact <- function(alignedFile, outputPath, referenceFile, incrArtefactThrEvery=10000000, verbosity=0) 
{
    
    ## Check that files to read exist
    
    # Aligned reads file
    if(!file.exists(alignedFile))
    {
        stop("The argument alignedFile does not refer to a valid file, please check that the path is correct and that the permissions policy allows to read it")
    }
    
    # Reference genome file (containing chromosomes names and lengths) 
    if(!file.exists(referenceFile))  
    {
        # If the file does not exist by itself, check if the user made reference to a precomputed one included in the package
        precomputedReferenceFilesFolder <- "resources" # This is the name of the folder (relative to the package installation) in which the precomputed reference files are stored

        # Try to get the full name of the reference file, returns empty string is not found 
        filePathToReferenceFile <- system.file(precomputedReferenceFilesFolder, paste(referenceFile, ".ref", sep="") ,package="Pasha")
        
        # If the file can't be found in the package, the user must have done something wrong
        if(0==nchar(filePathToReferenceFile))
        {
            # List the precomputed reference files available in the package
            precomputedReferenceGenomes <- gsub(".ref","",dir(system.file(precomputedReferenceFilesFolder ,package="Pasha"), pattern=".ref"))
            stop(paste("The argument referenceFile does not refer to a valid file or to a valid reference included in the package, if you provide your own file please check that the path is correct and that the permissions policy allows to read it. If you made reference to a file included in the package, please check spelling. Available ones :", paste(precomputedReferenceGenomes, collapse=" - ")))
        }
        else
        {
            # Let's use the precomputed reference file
            referenceFile <- filePathToReferenceFile
        }
    }
    
    
    ## Check other arguments consistency
        
    if(!(is.numeric(incrArtefactThrEvery) && (incrArtefactThrEvery>0)))
    {
        stop("The argument incrArtefactThrEvery must be a strictly positive numeric")
    }
    
    if(!(is.numeric(verbosity) && (verbosity>=0)))
    {
        stop("The argument verbosity must be a positive numeric")
    }

    
    ## Launch C function
    
    # Pre-allocate character string in memory for return
    sizeCharOutput <- 1024
    output <- paste(rep("_", sizeCharOutput),collapse="")
    
    # Arguments : (char** r_file_name, char** r_reference, int* r_incrArtefactThrEvery, int* r_verbosity)
    returnValues <- .C("C_RemoveArtifact", alignedFile, outputPath, referenceFile, as.integer(incrArtefactThrEvery), as.integer(verbosity), returnedOutputFile=output)

    return(returnValues$returnedOutputFile)
}
