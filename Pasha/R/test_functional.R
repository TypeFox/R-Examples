
# This function prepares the copy of the files in the directory where the test results will be written
# It then generate new results with the embeded SAM file
# Finally it compares this newly generated reult files with the precomputed ones that were previously unzipped
.testFunctional <- function(folderOutput=tempdir(), testType="regular", verbose=FALSE)
{
    if(!(file.exists(folderOutput) && file.info(folderOutput)$isdir)) stop("the path specified in folderOutput does not exist")
    if(!(testType %in% c("complete", "regular"))) stop("Invalid value for testType, must be \"complete\" or \"regular\"")
    
    message <- paste("Functional test of Pasha -", testType, "-")
    
    cat("\n\n")
    cat(rep("#",nchar(message)),"\n",sep="")
    cat(message,"\n")
    cat(rep("#",nchar(message)),"\n",sep="")
    
    generatedData_subfolderName <- "generated_Results"
    precomputedData_subfolderName <- "precomputed_Results"
    
    cat("\nPreparing folders and copying files for tests...")
    ### Create a subfolder for the test
    
    # generate a timestamp
    currentTime <- as.POSIXlt(Sys.time())
    stampPrefix <- format(currentTime, format="%y_%m_%d_%Hh%M")
    
    # create the subfolder where we will make two additional subfolders , one for computation of new results and one for the unzipping of precomputed results
    test_folderOutputGeneral <- file.path(folderOutput, paste("PashaTest_",stampPrefix,sep=""))
    dir.create(test_folderOutputGeneral, recursive=TRUE)
    
    test_folderOutputGenerated <- file.path(test_folderOutputGeneral, generatedData_subfolderName)
    dir.create(test_folderOutputGenerated, recursive=TRUE)
    
    test_folderOutputPrecomputed <- file.path(test_folderOutputGeneral, precomputedData_subfolderName)
    dir.create(test_folderOutputPrecomputed, recursive=TRUE)
    
    ### Copy original test data (BAM) to the folder and unzip the precomputed results
    
    # Get the path to and copy the embedded BAM file and its index
    testFileBAM_fileName <- "embedDataTest.bam"
    testFileBAM_fullPath <- system.file("extdata", testFileBAM_fileName,package="Pasha")

    testFileBAI_fileName <- "embedDataTest.bam.bai"
    testFileBAI_fullPath <- system.file("extdata", testFileBAI_fileName,package="Pasha")
    
    file.copy(testFileBAM_fullPath, test_folderOutputGenerated)
    file.copy(testFileBAI_fullPath, test_folderOutputGenerated)
    
    # Get the path to and copy the embedded multiloc data file
    testFileMulti_fileName <- "embededDataTest_MultiSignal.txt"
    testFileMulti_fullPath <- system.file("extdata", testFileMulti_fileName,package="Pasha")
    
    file.copy(testFileMulti_fullPath, test_folderOutputGenerated)
    
    
    ### unzip the precomputed expected results to the "precomputed" folder
    precomputedResultsZIP_fileName <- "resultsFromFunctionalTestsWIG.tar.gz"
    precomputedResultsZIP_fullPath <- system.file("extdata", precomputedResultsZIP_fileName,package="Pasha")
    
    untar(tarfile=precomputedResultsZIP_fullPath, compressed="gzip", exdir=test_folderOutputPrecomputed)# automatically create the subfolder precomputedResults
    
    cat("Done !\n")
    
    
    ### Generate results
    
    cat("\nStarting the results generation (silent output, see log file for details), this step can be long...")
    
    if(!verbose)
    {
        # Redirect R output temporarily to the log file to avoid endless log of Pasha pipeline
        sink(file=file.path(test_folderOutputGenerated,"resultGeneration_LogFile.log"))
    }
    
    # Generate results
    .testFunctional_generateResults(testFileBAM_folderName=test_folderOutputGenerated, testFileBAM_fileName=testFileBAM_fileName, testFileMulti_fileName=testFileMulti_fileName, testType=testType)
    
    if(!verbose)
    {
        # Stop the redirection of R output
        sink(NULL)
    }
    cat("Done !\n")
    
    ### Analyse differences
    cat("\nComparing the resulting files with precomputed ones...")
    
    # Get the filenames of all wigs and gff generated
    generatedResults_Files <- list.files(path=test_folderOutputGenerated, pattern="\\.wig$|\\.gff$|\\.bed$|\\.bw$", full.names=TRUE, recursive=TRUE)
    
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
.testFunctional_generateResults <- function(testFileBAM_folderName, testFileBAM_fileName, testFileMulti_fileName, testType="complete")
{
    
    if(!(testType %in% c("complete", "regular"))) stop("Invalid value for testType, must be \"complete\" or \"regular\"")
    
    # Prepare the path to the multiloc data (thet is copied in the same folder as the BAM)
    testFileMulti_fullPath <- file.path(testFileBAM_folderName, testFileMulti_fileName)
    
    
    # INPUT PARAMETERS
    INPUTFilesList <- list()
    INPUTFilesMultiList <- list()
    
    INPUTFilesList[["testBAM"]] <- list(folderName=testFileBAM_folderName, fileName=testFileBAM_fileName, fileType="BAM", chrPrefix="chr", chrSuffix="", pairedEnds=TRUE, midPoint=FALSE)
    INPUTFilesList[["testBAM_MIDPT_asSE"]] <- list(folderName=testFileBAM_folderName, fileName=testFileBAM_fileName, fileType="BAM", chrPrefix="chr", chrSuffix="", pairedEnds=FALSE, midPoint=TRUE)
    
    INPUTFilesMultiList[["testBAM"]] <- testFileMulti_fullPath
    INPUTFilesMultiList[["testBAM_MIDPT_asSE"]] <- testFileMulti_fullPath
    
    
    if(testType=="complete")
    {
        
        INPUTFilesList[["TestBAM_MIDPT"]] <- list(folderName=testFileBAM_folderName, fileName=testFileBAM_fileName, fileType="BAM", chrPrefix="chr", chrSuffix="", pairedEnds=TRUE, midPoint=TRUE)
        INPUTFilesList[["TestBAM_asSE"]] <- list(folderName=testFileBAM_folderName, fileName=testFileBAM_fileName, fileType="BAM", chrPrefix="chr", chrSuffix="", pairedEnds=FALSE, midPoint=FALSE)
        
        INPUTFilesMultiList[["TestBAM_MIDPT"]] <- testFileMulti_fullPath
        INPUTFilesMultiList[["TestBAM_asSE"]] <- testFileMulti_fullPath
    }
    
    
    suppressWarnings(
    processPipeline(INPUTFilesList=INPUTFilesList,
            resultSubFolder              = "Results",
            reportFilesSubFolder         = "Reports",
            WIGfs                        = TRUE,
            WIGvs                        = TRUE,
            GFF                          = FALSE, # GFF is not tested since the resulting files are too big to be stored in package
            BED                          = TRUE,
			BIGWIG                       = FALSE, # BIGWIG is not tested since the resulting files are too big to be stored in package
            compatibilityOutputWIG       = FALSE,
            # COMPLEX PARAMETERS (ATOMIC OR VECTORS OR LIST OF IT)
            incrArtefactThrEvery         = if(testType=="complete") c(-2,NA)      else c(-2),
            binSize                      = if(testType=="complete") c(1, 50, 200) else c(200),
            elongationSize               = if(testType=="complete") c(NA, 100, 0) else c(NA),
            rangeSelection               = if(testType=="complete") c(2)          else IRanges(0,-1),
            annotationFilesGFF           = NA,
            annotationGenomeFiles        = NA,
            # ATOMIC PARAMETERS
            elongationEstimationRange    = c(mini=30, maxi=400, by=10),
            rehabilitationStep           = c("orphans","orphansFromArtefacts"),
            removeChrNamesContaining     = "random|hap",
            ignoreInsertsOver            = 500,
            nbCPUs                       = 1,
            keepTemp                     = TRUE, # Keep the intermediary files that led to the final ones (rehab and multi)
            logTofile                    = NULL,
            eraseLog                     = TRUE,
            # LIST PARAMETERS (one element per expName)
            multiLocFilesList            = INPUTFilesMultiList)) # A list with experiments names and associated filenames to treat
    
    
}