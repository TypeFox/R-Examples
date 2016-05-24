##' zcrana function adapted from libassp
##'
##' Analysis of the averages of the short-term positive and
##' negative zero-crossing rates of the signal in <listOfFiles>.
##' Analysis results will be written to a file with the
##' base name of the input file and extension '.zcr'.
##' Default output is in SSFF binary format (track 'zcr').
##' @title zcrana
##' @param listOfFiles vector of file paths to be processed by function 
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default: begin of file)
##' @param centerTime = <time>  set single-frame analysis with the analysis window centred at <time> seconds; 
##' overrules beginTime, endTime and windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default: end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize = <dur>:  set analysis window size to <dur> ms (default: 25.0)
##' @param toFile write results to file (default extension is .zcr)
##' @param explicitExt set if you wish to overwride the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param forceToLog is set by the global package variable useWrasspLogger. This is set
##' to FALSE by default and should be set to TRUE is logging is desired.
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @useDynLib wrassp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate zcr values
##' res <- zcrana(path2wav, toFile=FALSE)
##' 
##' # plot zcr values
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res$zcr, 
##'      type='l', 
##'      xlab='time (s)', 
##'      ylab='ZCR values')
##' 
##' @export
'zcrana' <- function(listOfFiles = NULL, optLogFilePath = NULL, 
                     beginTime = 0.0, centerTime = FALSE, 
                     endTime = 0.0, windowShift = 5.0, 
                     windowSize = 25.0, toFile = TRUE, 
                     explicitExt = NULL, outputDirectory = NULL,
                     forceToLog = useWrasspLogger){
  
  ###########################
  # a few parameter checks and expand paths
  
  if (is.null(listOfFiles)) {
    stop(paste("listOfFiles is NULL! It has to be a string or vector of file",
               "paths (min length = 1) pointing to valid file(s) to perform",
               "the given analysis function."))
  }
  
  if (is.null(optLogFilePath) && forceToLog){
    stop("optLogFilePath is NULL! -> not logging!")
  }else{
    if(forceToLog){
      optLogFilePath = path.expand(optLogFilePath)  
    }
  }
  
  if (!is.null(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory))
    finfo  <- file.info(outputDirectory)
    if (is.na(finfo$isdir))
      if (!dir.create(outputDirectory, recursive=TRUE))
        stop('Unable to create output directory.')
    else if (!finfo$isdir)
      stop(paste(outputDirectory, 'exists but is not a directory.'))
  }
  ###########################
  # Pre-process file list
  listOfFiles <- prepareFiles(listOfFiles)
  
  ###########################
  # perform analysis
  
  if(length(listOfFiles)==1){
    pb <- NULL
  }else{
    if(toFile==FALSE){
      stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
    }
    cat('\n  INFO: applying zcrana to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }	
  
  externalRes = invisible(.External("performAssp", PACKAGE = "wrassp", 
                                    listOfFiles, fname = "zcrana", 
                                    beginTime = beginTime, centerTime = centerTime, 
                                    endTime = endTime, windowShift = windowShift, 
                                    windowSize = windowSize, 
                                    toFile = toFile, explicitExt = explicitExt, 
                                    outputDirectory = outputDirectory, progressBar = pb))
  
  
  ############################
  # write options to options log file
  if (forceToLog){
    optionsGivenAsArgs = as.list(match.call(expand.dots = TRUE))
    wrassp.logger(optionsGivenAsArgs[[1]], optionsGivenAsArgs[-1],
                  optLogFilePath, listOfFiles)
  }    
  
  #############################
  # return dataObj if length only one file
  
  if(!(length(listOfFiles)==1)){
    close(pb)
  }else{
    return(externalRes)
  }
}

