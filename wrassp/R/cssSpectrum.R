##' calculate cepstrally smoothed spectrum using libassp
##'
##' Short-term spectral analysis of the signal in <listOfFiles>
##' using the Fast Fourier Transform and cepstral smoothing.
##' Analysis results will be written to a file with the
##' base name of the input file and '.css.' as extension.
##' Default output is in SSFF format with
##' 'css' in lower case as track name.
##' @title cssSpectrum
##' @param listOfFiles vector of file paths to be processed by function 
##' @param optLogFilePath path to option log file
##' @param beginTime = <time>: set begin of analysis interval to <time> seconds
##' (default: begin of data)
##' @param centerTime = <time>: set single-frame analysis with the analysis
##' window centred at <time> seconds; overrules beginTime, endTime and
##' windowShift options
##' @param endTime = <time>: set end of analysis interval to <time> seconds
##' (default: end of data)
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a freqequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##' @param windowShift = <dur>: set analysis window shift to <dur> ms
##' (default: 5.0)
##' @param window = <type>: set analysis window function to <type> (default:
##' BLACKMAN)
##' @param numCeps = <num>: set number of cepstral coeffcients used to <num>
##' (default: sampling rate in kHz + 1; minimum: 2)
##' @param toFile write results to file (default extension depends on )
##' @param explicitExt set if you wish to overwride the default extension
##' @param outputDirectory directory in which output files are stored. Defaults to NULL, i.e.
##' the directory of the input files
##' @param forceToLog is set by the global package variable useWrasspLogger. This is set
##' to FALSE by default and should be set to TRUE is logging is desired.
##' @return nrOfProcessedFiles or if only one file to process return
##' AsspDataObj of that file
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @seealso \code{\link{dftSpectrum}}, \code{\link{lpsSpectrum}}, \code{\link{cepstrum}}; 
##' all derived from libassp's spectrum function.
##' @useDynLib wrassp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"), 
##'                        pattern = glob2rx("*.wav"), 
##'                        full.names = TRUE)[1]
##' 
##' # calculate cepstrally smoothed spectrum
##' res <- cssSpectrum(path2wav, toFile=FALSE)
##' 
##' # plot spectral values at midpoint of signal
##' plot(res$css[dim(res$css)[1]/2,], 
##'      type='l', 
##'      xlab='spectral value index', 
##'      ylab='spectral value')
##'      
##' @export
'cssSpectrum' <- function(listOfFiles = NULL, optLogFilePath = NULL,
                          beginTime = 0.0, centerTime = FALSE,
                          endTime = 0.0, resolution = 40.0,
                          fftLength = 0, windowShift = 5.0, 
                          window = 'BLACKMAN', numCeps = 0, 
                          toFile = TRUE, explicitExt = NULL, 
                          outputDirectory = NULL, forceToLog = useWrasspLogger){
  
  ## ########################
  ## a few parameter checks and expand paths
  
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
  
  ## #######################
  ## perform analysis
  
  if(length(listOfFiles)==1){
    pb <- NULL
  }else{
    if(toFile==FALSE){
      stop("length(listOfFiles) is > 1 and ToFile=FALSE! ToFile=FALSE only permitted for single files.")
    }
    cat('\n  INFO: applying cssSpectrum to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }	
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "spectrum", beginTime = beginTime, 
                                    centerTime = centerTime, endTime = endTime, 
                                    spectrumType = 'CSS',
                                    resolution = resolution, 
                                    fftLength = as.integer(fftLength), 
                                    windowShift = windowShift, window = window, 
                                    numCeps = as.integer(numCeps), 
                                    toFile = toFile, explicitExt = explicitExt, 
                                    progressBar = pb, outputDirectory = outputDirectory,
                                    PACKAGE = "wrassp"))
  
  
  ## #########################
  ## write options to options log file
  if (forceToLog){
    optionsGivenAsArgs = as.list(match.call(expand.dots = TRUE))
    wrassp.logger(optionsGivenAsArgs[[1]], optionsGivenAsArgs[-1],
                  optLogFilePath, listOfFiles)
  }
  
  ## #########################
  ## return dataObj if length only one file
  
  if(!(length(listOfFiles)==1)){
    close(pb)
  }else{
    return(externalRes)
  }
}
