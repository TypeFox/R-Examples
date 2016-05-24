##' acfana function adapted from libassp
##'
##' Analysis of short-term autocorrelation function of
##' the signals in <listOFFiles>.
##' Analysis results will be written to a file with the
##' base name of the input file and extension '.acf'.
##' Default output is in SSFF binary format (track 'acf').
##' @title acfana
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default: 0 = beginning of file)
##' @param centerTime = <time>: set single-frame analysis with the analysis window centred at <time> seconds; 
##' overrules BeginTime, EndTime and WindowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default: 0 = end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize = <dur>: set analysis window size to <dur> ms; overrules EffectiveLength parameter
##' @param effectiveLength make window size effective rather than exact
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param analysisOrder = <num>: set analysis order to <num> (default: 0 = sample rate in kHz + 3)
##' @param energyNormalization calculate energy-normalized autocorrelation
##' @param lengthNormalization calculate length-normalized autocorrelation
##' @param toFile write results to file (default extension is .acf)
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
##' # calculate short-term autocorrelation
##' res <- acfana(path2wav, toFile=FALSE)
##' 
##' # plot short-term autocorrelation values
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) + 
##'         attr(res, 'startTime'), 
##'         res$acf, 
##'         type='l', 
##'         xlab='time (s)', 
##'         ylab='short-term autocorrelation values')
##'         
##' @export
'acfana' <- function(listOfFiles = NULL, optLogFilePath = NULL, 
                     beginTime = 0.0, centerTime = FALSE, 
                     endTime = 0.0, windowShift = 5.0, 
                     windowSize = 20.0, effectiveLength = TRUE, 
                     window = "BLACKMAN", analysisOrder = 0, 
                     energyNormalization = FALSE, lengthNormalization = FALSE, 
                     toFile = TRUE, explicitExt = NULL, outputDirectory = NULL,
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
  
  if(!isAsspWindowType(window)){
    stop("WindowFunction of type '", window,"' is not supported!")
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
    cat('\n  INFO: applying acfana to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "acfana", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    windowShift = windowShift, windowSize = windowSize, 
                                    effectiveLength = effectiveLength, window = window, 
                                    analysisOrder = as.integer(analysisOrder), energyNormalization = energyNormalization, 
                                    lengthNormalization = lengthNormalization, toFile = toFile, 
                                    explicitExt = explicitExt, progressBar = pb,
                                    outputDirectory = outputDirectory, PACKAGE = "wrassp"))
  
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
