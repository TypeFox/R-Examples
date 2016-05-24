##' afdiff function adapted from libassp
##'
##' Computes the first difference of the signal in the audio-
##' formatted file(s) <listOfFiles>. The differentiated signal will
##' be written to a file with the base name of the input file
##' and an extension consisting of '.d', followed by the
##' extension of the input file. The format of the output file
##' will be the same as that of the input file.
##' Differentiation can improve results on F0 analysis of e.g.
##' EGG signals because it removes a DC offset, attenuates
##' very low frequency components - and in the case of central
##' differentiation also very high ones - and enhances the
##' moment of glottal closure.
##' @title afdiff
##' @param listOfFiles vector of file paths to be processed by function
##' @param optLogFilePath path to option log file
##' @param computeBackwardDifference compute backward difference (s'[n] = s[n] - s[n-1]) (default: forward difference s'[n] = s[n+1] - s[n])
##' @param computeCentralDifference compute central/interpolated/3-point difference
##' @param channel = <num>: for multi-channel input files: extract and differentiate channel <num> (1 <= <num> <= 8  default: channel 1)
##' @param toFile write results to file (default extension is .d+(extensionsOfAudioFile))
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
##' # compute the first forward difference of the signal
##' res <- afdiff(path2wav, toFile=FALSE)
##' 
##' # plot samples
##' # (only plot every 10th element to accelerate plotting)
##' plot(seq(0,numRecs.AsspDataObj(res) - 1, 10) / rate.AsspDataObj(res), 
##'      res$audio[c(TRUE, rep(FALSE,9))], 
##'      type='l', 
##'      xlab='time (s)', 
##'      ylab='Audio samples')
##'      
##' @export
'afdiff' <- function(listOfFiles = NULL, optLogFilePath = NULL, 
                     computeBackwardDifference = FALSE, computeCentralDifference = FALSE, 
                     channel = 1, toFile = TRUE, 
                     explicitExt=NULL, outputDirectory = NULL,
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
    cat('\n  INFO: applying afdiff to', length(listOfFiles), 'files\n')
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }
  
  externalRes = invisible(.External("performAssp", listOfFiles, 
                                    fname = "afdiff", computeBackwardDifference = computeBackwardDifference, 
                                    channel = as.integer(channel), toFile = toFile, 
                                    explicitExt = explicitExt, progressBar=pb, 
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
