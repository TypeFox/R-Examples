################################################################################
#' Analyze power band by episodes
#' @description
#' Analyzes the ULF, VLF, LF and HF bands from a given indexFreqAnalysis allowing
#'  to evaluate the application of a desired function inside and outside each episode.
#' @param HRVData Data structure that stores the beats register and information related to it. 
#' @param indexFreqAnalysis Integer value denoting which frequency analysis is going to be analyzed using func. Default: 1
#' @param Tag Type of episode
#' @param verbose Deprecated argument maintained for compatibility, use SetVerbose() instead
#' @param func Function to be applied to each power band inside and outside episodes
#' @param ... Optional arguments for func.
#' @return Returns a list with two objects, that is, the values of the application of the selected function
#'  inside ("resultIn") and outside ("resultOut") episodes in the given indexFreqAnalysis. Each of these 
#'  list has another set of lists: the  "ULF", "VLF", "LF" and "HF" lists.
#' @examples 
#' \dontrun{
#' hrv.data = CreateHRVData()
#' hrv.data = SetVerbose(hrv.data, TRUE)
#' hrv.data = LoadBeat(hrv.data, fileType = "WFDB", "a03", RecordPath ="beatsFolder/", 
#'                     annotator = "qrs")
#'                     hrv.data = LoadApneaWFDB(hrv.data, RecordName="a03",Tag="Apnea",
#'                                              RecordPath="beatsFolder/")
#' hrv.data = BuildNIHR(hrv.data)
#' hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
#' hrv.data = CreateFreqAnalysis(hrv.data)
#' hrv.data = CalculatePowerBand( hrv.data , indexFreqAnalysis= 1,
#'                                type = "wavelet", wavelet = "la8",
#'                                 bandtolerance = 0.01, relative = FALSE)
#' results = AnalyzePowerBandsByEpisodes(hrv.data,indexFreqAnalysis=1,
#'                                        Tag="Apnea",func=mean)}
AnalyzePowerBandsByEpisodes = function(HRVData, indexFreqAnalysis = length(HRVData$FreqAnalysis), Tag="", verbose=NULL,func, ...) {
    # ----------------------------------------------
    # Analyzes PowerBands using Episodes information
    # ----------------------------------------------
    #  indexFreqAnalysis -> which frequency analysis is going to be analyzed using func
    #  Tag -> specifies tag of episodes
    #  func -> function to apply 
    #  ... -> additional arguments for func
    #  Returns a list with two objects result
    
    funcToApply =  match.fun(func)
    nameFunc = deparse(substitute(func))
    
    
    #check if indexFreqAnalysis exists
    if ((length(HRVData$FreqAnalysis) < indexFreqAnalysis) || (indexFreqAnalysis<1) ) {
      stop("   --- Frequency analysis no. ",indexFreqAnalysis," not present!! ---\n    --- Quitting now!! ---\n")
    }
    
    if (!is.null(verbose)) {
      cat("  --- Warning: deprecated argument, using SetVerbose() instead ---\n    --- See help for more information!! ---\n")
      SetVerbose(HRVData,verbose)
    }
    
    if (HRVData$Verbose) {
      cat("** Applying function to power bands in frequency analysis" ,indexFreqAnalysis," using episodic information **\n");
      cat("   Function: ",nameFunc,"()\n",sep="")
    }
    
    if (is.null(HRVData$Episodes)) {
      stop("  --- Episodes not present\n    --- Quitting now!! ---\n")
    }
    
    if (HRVData$Verbose) {
      if (Tag=="") {
        cat("   No tag was specified\n")
      } else {
        cat("   Using episodes with tag:",Tag,"\n")
      }
    }
    
    episodicInformation = SplitPowerBandByEpisodes(HRVData,
                                                   indexFreqAnalysis = indexFreqAnalysis,
                                                   Tag = Tag)
    
    bandNames = names(episodicInformation$InEpisodes)
    resultIn = list()
    resultOut = list()
    for (band in bandNames){
      resultIn[[band]] = funcToApply(episodicInformation$InEpisodes[[band]], ...)
      resultOut[[band]] = funcToApply(episodicInformation$OutEpisodes[[band]], ...)
    }
  
    
    result=list(resultIn=resultIn,resultOut=resultOut)
    return(result)

}
