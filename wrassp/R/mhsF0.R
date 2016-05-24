##' mhsF0 function adapted from libassp
##'
##' Pitch analysis of the speech signal in <listOfFile> using
##' Michel's/Modified Harmonic Sieve algorithm.
##' Analysis results will be written to a file with the
##' base name of the input file and extension '.pit'.
##' Default output is in SSFF binary format (track 'pitch').
##' @title mhsF0
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds (default = 0: begin of file) 
##' @param centerTime = <time>:  set single-frame analysis with the analysis 
##' window centred at <time> seconds; overrules beginTime, endTime and windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds (default = 0: end of file)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms (default: 5.0)
##' @param gender = <code>  set gender-specific pitch ranges; <code> may be:
##' "f[emale]" (80.0 - 600.0 Hz)
##' "m[ale]" (50.0 - 375.0 Hz)
##' "u[nknown]" (default; 50.0 - 600.0 Hz)
##' @param maxF = <freq>: set maximum pitch value to <freq> Hz (default: 500.0)
##' @param minF = <freq>:  set minimum pitch value to <freq> Hz (default: 50.0  minimum: 25.0)
##' @param minAmp = <amp>:  minimum signal amplitude (default: 50)
##' @param minAC1 = <freq>: minimum 1st correlation coefficient (default: 0.250)
##' @param minRMS = <num>:  minimum RMS amplitude in dB (default: 18.0)
##' @param maxZCR = <freq>: maximum zero crossing rate in Hz (default: 3000)
##' @param minProb = <num>: minimum quality value of F0 fit (default: 0.520)
##' @param plainSpectrum use plain rather than masked power spectrum
##' @param toFile write results to file (default extension is .pit)
##' @param explicitExt set if you wish to overwride the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param forceToLog is set by the global package variable useWrasspLogger. This is set
##' to FALSE by default and should be set to TRUE is logging is desired.
##' @return nrOfProcessedFiles or if only one file to process return AsspDataObj of that file
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @aliases mhspitch f0_mhs
##' @seealso \code{\link{ksvF0}} for an tracking the fundamental frequency
##' @useDynLib wrassp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate fundamental frequency contour
##' res <- mhsF0(path2wav, toFile=FALSE)
##' 
##' # plot fundamental frequency contour
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res$pitch, 
##'      type='l', 
##'      xlab='time (s)', 
##'      ylab='F0 frequency (Hz)')
##' 
##' @export
'mhsF0' <- 'mhspitch' <- 'f0_mhs' <-function(listOfFiles = NULL, optLogFilePath = NULL,
                                             beginTime = 0.0, centerTime = FALSE, 
                                             endTime = 0.0, windowShift = 5.0, 
                                             gender = 'u', maxF = 600.0, 
                                             minF = 50.0, minAmp = 50.0, 
                                             minAC1 = 0.25, minRMS = 18.0, 
                                             maxZCR = 3000.0, minProb = 0.52, 
                                             plainSpectrum = FALSE, toFile = TRUE, 
                                             explicitExt = NULL,  outputDirectory = NULL,
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
    OutputDirectory = normalizePath(path.expand(outputDirectory))
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
    cat('\n  INFO: applying mhspitch to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }		
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "mhspitch", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    windowShift = windowShift, gender = gender, 
                                    maxF = maxF, minF = minF, 
                                    minAmp = minAmp, minAC1 = minAC1, 
                                    minRMS = minRMS, maxZCR = maxZCR, 
                                    minProb = minProb, plainSpectrum = plainSpectrum, 
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
