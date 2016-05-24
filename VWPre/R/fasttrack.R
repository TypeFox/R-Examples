#' Fast-track processing
#' 
#' \code{fasttrack} is a meta-function for advanced users who are 
#' already familiar with the package functions.
#' It takes all necessary arguments for the component functions
#' to produce proportion looks. It also takes an argument specifying 
#' if empirical logits or binomial data should be output.  The
#' function returns a dataframe contained the result of the series 
#' of subroutines.
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' 
#' @param data A data frame object created from an Eyelink Sample Report.
#' @param Subject An obligatory string containing the column name corresponding 
#' to the subject identifier.
#' @param Item An optional string containing the column name corresponding to 
#' the item identifier; by default, NA.
#' @param EventColumns A vector specifying the columns which will be used for creating 
#' the Event variable; by default, Subject and TRIAL_INDEX.
#' @param NoIA A positive integer indicating the number of interest areas defined 
#' when creating the study. 
#' @param Offset A positive integer indicating amount of time in milliseconds.
#' @param Recording A string indicating which eyes were used for recording gaze data.
#' @param WhenLandR A string indicating which eye ("Right" or "Left) to use 
#' if gaze data is available for both eyes (i.e., Recording = "LandR"). 
#' @param BinSize A positive integer indicating the size of the binning window 
#' (in milliseconds).
#' @param SamplingRate A positive integer indicating the sampling rate (in Hertz) 
#' used to record the gaze data.
#' @param SamplesPerBin A positive integer indicating the number of samples in
#' each bin.
#' @param Constant A positive number used for the empirical logit and weights
#' calculation; by default, 0.5 as in Barr (2008).
#' @param Output An obligatory string containing either "ELogit" or "Binomial".
#' @return A data table containing formatting and calculations.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Perform meta-function on data
#' df <- fasttrack(data = dat, Subject = "RECORDING_SESSION_LABEL", Item = "itemid", 
#'        EventColumns = c("Subject", "TRIAL_INDEX"), NoIA = 4, Offset = 100, 
#'				Recording = "LandR", WhenLandR = "Right", BinSize = 20, 
#'				SamplingRate = 1000, SamplesPerBin = 20, Constant = 0.5, 
#'				Output = "ELogit")
#' }
fasttrack = function(data = data, Subject = Subject, Item = Item, 
                        EventColumns = c("Subject", "TRIAL_INDEX"), NoIA = NoIA,
                        Offset = Offset, Recording = Recording, 
                        WhenLandR = WhenLandR, BinSize = BinSize, SamplingRate = SamplingRate,
                        SamplesPerBin = SamplesPerBin, Constant = 0.5, Output = Output) {
  
  dat <- data
  Subject <- Subject
  Item <- Item
  EventColumns <- EventColumns
  NoIA <- NoIA
  Offset <- Offset
  Recording <- Recording
  WhenLandR <- WhenLandR
  BinSize <- BinSize
  SamplingRate <- SamplingRate
  SamplesPerBin <- SamplesPerBin
  Constant <- Constant
  Output <- Output
  
  print("Preparing data...")
  dat0 <- prep_data(data = dat, Subject = Subject, Item = Item, EventColumns = EventColumns)

  print(paste("Relabelling outside of", NoIA, "interest areas...", sep = " "))
  dat1 <- relabel_na(data = dat0, NoIA = NoIA)
  rm(dat0)
  
  print(paste("Creating time series with", Offset, "ms offset...", sep = " "))
  dat2 <- create_time_series(data = dat1, Offset = Offset)
  rm(dat1)
  
  check_time_series(data = dat2)
  check_eye_recording(data = dat2)
  
  print(paste("Selecting recorded eye..."))
  dat3 <- select_recorded_eye(data = dat2, Recording = Recording, WhenLandR = WhenLandR)
  rm(dat2)
  
  check_samplingrate(dat3)
  
  print(paste("Binning", SamplingRate, "Hz data into", BinSize, "ms bins..."))
  print("Calculating proportions...")
  dat4 <- bin_prop(dat3, NoIA = NoIA, BinSize = BinSize, SamplingRate = SamplingRate)
  rm(dat3)
  
  check_samplingrate(dat4)
  check_samples_per_bin(dat4)
  
  print(paste("Preparing", Output, "output...", sep = " "))
  if (Output == "ELogit") {
    dat5 <- transform_to_elogit(dat4, NoIA = 4, SamplesPerBin = 20, Constant = Constant)
  } else if (Output == "Binomial") {
    dat5 <- create_binomial(data = dat4, NoIA = 4)
  }
  
  rm(dat4)
  
  return(dat5)
  
}