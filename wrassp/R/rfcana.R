##' rfcana function adapted from libassp 
##'
##' Linear Prediction analysis of <listOfFiles> using the
##' autocorrelation method and the Durbin recursion.
##' This program calculates the RMS amplitudes of the input
##' and residual signal in dB and, per default, reflection
##' coefficients (see '-t' option).
##' Analysis results will be written to a file with the
##' base name of the input file and the parameter type in
##' lower case as extension (e.g. '.rfc').
##' Default output is in SSFF binary format (tracks 'rms',
##' 'gain' and the LP type in lower case).
##' @title rfcana
##' @param listOfFiles vector of file paths to be processed by function 
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default = 0: begin of file)
##' @param centerTime set single-frame analysis with the analysis window centred at <time> seconds; 
##' overrules beginTime, endTime and windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default = 0: end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param windowSize = <dur>: set analysis window size to <dur> ms; overrules effectiveLength option
##' @param effectiveLength make window size effective rather than exact
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param order = <num>: set prediction order to <num> (default: sample rate in kHz + 3)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (default: -0.95)
##' @param lpType = <type>: calculate <type> LP parameters; <type> may be:
##' "ARF": area function
##' "LAR": log area ratios
##' "LPC": linear prediction filter coefficients
##' "RFC": reflection coefficients (default)
##' @param toFile write results to file (default extension dependent on LpType .arf/.lar/.lpc/.rfc)
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
##' # perform linear prediction analysis
##' res <- rfcana(path2wav, toFile=FALSE)
##' 
##' # plot reflection coefficients
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) + 
##'           attr(res, 'startTime'), 
##'         res$rfc, 
##'         type='l', 
##'         xlab='time (s)', 
##'         ylab='reflection coefficient values')
##'         
##' @export
'rfcana' <- function(listOfFiles = NULL, optLogFilePath = NULL, 
                     beginTime = 0.0, centerTime = FALSE, 
                     endTime = 0.0, windowShift = 5.0, 
                     windowSize = 20.0, effectiveLength = TRUE, 
                     window = 'BLACKMAN', order = 0, 
                     preemphasis = -0.95, lpType = 'RFC', 
                     toFile = TRUE, explicitExt = NULL,
                     outputDirectory = NULL, forceToLog = useWrasspLogger){
  
  
  ###########################
  # a few parameter checks and expand files
  
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
  
  if(!isAsspLpType(lpType)){
    stop("LpType of type '", lpType,"' is not supported!")
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
    cat('\n  INFO: applying rfcana to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }
  
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "rfcana", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    windowShift = windowShift, windowSize = windowSize, 
                                    effectiveLength = effectiveLength, window = window, 
                                    order = as.integer(order), 
                                    preemphasis = preemphasis, lpType = lpType, 
                                    toFile = toFile, explicitExt = explicitExt, 
                                    progressBar = pb, outputDirectory = outputDirectory,
                                    PACKAGE = "wrassp"))
  
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
