
# This function prepares the copy of the files in the directory where the test results will be written
# It then generate new results with the embeded bow file
# Finally it compares this newly generated result files with the precomputed ones that were previously unzipped
.testFunctionalMultiread <- function(folderOutput=tempdir(), testType="regular", verbose=FALSE)
{
    if(!(file.exists(folderOutput) && file.info(folderOutput)$isdir)) stop("the path specified in folderOutput does not exist")
    if(!(testType %in% c("complete", "regular"))) stop("Invalid value for testType, must be \"complete\" or \"regular\"")
    
    message <- paste("Functional test of Pasha (Multiread) -", testType, "-")
    
    cat("\n\n")
    cat(rep("#",nchar(message)),"\n",sep="")
    cat(message,"\n")
    cat(rep("#",nchar(message)),"\n",sep="")
    
    generatedData_subfolderName <- "generated_MultireadResults"
    precomputedData_subfolderName <- "precomputed_MultireadResults"
    
    cat("\nPreparing folders and copying files for tests...")
    ### Create a subfolder for the test
    
    # generate a timestamp
    currentTime <- as.POSIXlt(Sys.time())
    stampPrefix <- format(currentTime, format="%Y_%m_%d_%Hh%M")
    
    # create the subfolder where we will make two additional subfolders , one for computation of new results and one for the unzipping of precomputed results
    test_folderOutputGeneral <- file.path(folderOutput, paste("R_library_Pasha_multireadfunctional_test_",stampPrefix,sep=""))
    dir.create(test_folderOutputGeneral, recursive=TRUE)
    
    test_folderOutputGenerated <- file.path(test_folderOutputGeneral, generatedData_subfolderName)
    dir.create(test_folderOutputGenerated, recursive=TRUE)
    
    test_folderOutputPrecomputed <- file.path(test_folderOutputGeneral, precomputedData_subfolderName)
    dir.create(test_folderOutputPrecomputed, recursive=TRUE)
    
    ### Copy original test data (BOW) to the folder and unzip the precomputed results
    
    # Get the path to and copy the embedded BOW file
    testFileBOW_fileName <- "embededDataTest_MultiSignal.bow"
    testFileBOW_fullPath <- system.file("extdata", testFileBOW_fileName,package="Pasha")
    
    file.copy(testFileBOW_fullPath, test_folderOutputGenerated)
    
    ### unzip the precomputed expected results to the "precomputed" folder
    precomputedResultsZIP_fileName <- "resultsFromMultireadFunctionalTests.tar.gz"
    precomputedResultsZIP_fullPath <- system.file("extdata", precomputedResultsZIP_fileName,package="Pasha")
    
    untar(precomputedResultsZIP_fullPath, compressed="gzip", exdir=test_folderOutputPrecomputed)# automatically create the subfolder precomputedResults
    
    cat("Done !\n")
    
    ### Generate results
    
    cat("\nStarting the results generation (silent output, see log file for details), this step can be long...")
    
    verbosity  <-  1
    if(!verbose)
    {
        # Redirect R output temporarily to the log file to avoid endless log of Pasha pipeline
        sink(file=file.path(test_folderOutputGenerated,"resultMultireadGeneration_LogFile.log"))
        verbosity  <-  0
    }
    
    # Generate results
    .testFunctionalMultiread_generateResults(testFileBOW_folderName=test_folderOutputGenerated, testFileBOW_fileName=testFileBOW_fileName, folderOutputGenerated = test_folderOutputGenerated, verbosity)
    
    if(!verbose)
    {
        # Stop the redirection of R output
        sink(NULL)
    }
    cat("Done !\n")
    
    ### Analyse differences
    cat("\nComparing the resulting files with precomputed ones...")
    
    # Get the filenames of all files we want to compare to reference
    generatedResults_Files <- list.files(path=test_folderOutputGenerated, pattern="txt|bow", full.names=TRUE, recursive=TRUE)

    # Remove the ones that we want to exclude from regular tests (these test are started automatically and might prevent the package validation if they fail. This is likely if the float representation is different on different machines/platforms, they should be started manually on demand only)
    if(testType!="complete") generatedResults_Files <- generatedResults_Files[!grepl("*Dispatch*", generatedResults_Files)]

    cat("\n\nNumber of results files generated :",length(generatedResults_Files),"\n")
    
    # Generate corresponding names for precomputed subfolder
    precomputedResults_Files <- gsub(generatedData_subfolderName,precomputedData_subfolderName,generatedResults_Files)
    
    if(!all(file.exists(precomputedResults_Files))) stop("Some files generated can't be found in the precomputed file list, make sure that the precomputed file archive is up to date or that the parameters of tests are compatible.")
    
    # Compute the md5sum of the files and compare them
    md5_generatedResults <- md5sum(generatedResults_Files)
    md5_precomputedResults <- md5sum(precomputedResults_Files)
    
    comparison_md5_OK <- (md5_generatedResults==md5_precomputedResults)
    
    cat("\nNumber of successful comparison (md5sum) with precomputed files : ",sum(comparison_md5_OK),"/", length(comparison_md5_OK), sep="")
    
    if(!all(comparison_md5_OK))
    {
        cat("\n\nDifferences were found between at least one generated file and its precomputed version :\n")
        cat(paste(names(md5_generatedResults[!comparison_md5_OK]),collapse="\n"))
        
        stop("Error while comparing results generated and precomputed bundled results")
    }
    
    cat("\n\nDone, all results seem consistent with reference\n\n")
    
}

# This function gnerate the results for a test "complete" or "regular" that aim to go through all the combination of parameters
.testFunctionalMultiread_generateResults <- function(testFileBOW_folderName, testFileBOW_fileName, folderOutputGenerated, verbosity=0)
{
    # Define the input file
    input_file_BOW  <-  file.path( folderOutputGenerated, testFileBOW_fileName)
    
    # Execute the Dispatching tests
    cat("\n\nAllocating scores using Uniform method")
    multiread_UniformDispatch(input_file_BOW, folderOutputGenerated, "mm9", 8, verbosity)
    
    cat("\n\nAllocating scores using CSEM method")
    multiread_CSEMDispatch(input_file_BOW, folderOutputGenerated, "mm9", 101, 200, 8, verbosity)
}