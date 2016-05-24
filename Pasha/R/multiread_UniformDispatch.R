
multiread_UniformDispatch <- function(alignedFile, outputFolder, referenceFile, incrArtefactThrEvery = NA, verbosity = 0) 
{
    
    # Check if the alignedFile exists
    if(!file.exists(alignedFile)){
        stop("The argument alignedFile does not refer to a valid file, please check that the path is correct and that the permissions policy allows to read it")
    }
    
    # Check the existence of output folder. If it exists and is not a valid dir, stop script.
    # If it does not exists, create it
    .safeCreateFolder( outputFolder)
    
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
    
    
    ## Launch C function
    sizeCharOutput <- 1024
    
    # Remove artifact if asked
    if( !is.na( incrArtefactThrEvery)){
        if( is.numeric(incrArtefactThrEvery) & incrArtefactThrEvery > 0){
            # Execute the artifact removal
            alignedFile  <-  multiread_RemoveArtifact(alignedFile, outputFolder, referenceFile, incrArtefactThrEvery, verbosity)
            # Check if the returned alignedFile exists
            if(!file.exists(alignedFile)){
                stop(paste("The artifact removal encountered an issue and does not return a correct output file:", alignedFile))
            }
        } else{
            stop("The argument 'incrArtefactThrEvery' must be a strictly positive number.")
        }    
    }
        
    # Arguments : (char** r_file_name, char** r_output_dir, char** r_reference)
    output <- paste(rep("_", sizeCharOutput),collapse="")
    returnValues <- .C("C_UniformDispatch", alignedFile, outputFolder, referenceFile, returnedOutputFile = output)
    
    return(returnValues$returnedOutputFile)
}
