#' Check the new time series
#' 
#' \code{check_time_series} examines the first value in the Time column
#' for each event. If they are equal, it will return a single value.  The returned
#' value should be equal to 0 minus the offset.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{create_time_series}}.
#' @return The value(s) of Time (in milliseconds) at which events begin relative 
#' to the onset of the auditory stimulus.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Check the starting Time column...
#' check_time_series(dat)
#' }
check_time_series = function(data = data) {
  event_start_table = data %>%
    summarise(ftime = min(Time))
  print(unique(event_start_table$ftime))
}


#' Check the number of samples in each bin
#' 
#' \code{check_samples_per_bin} determines the number of samples in each
#' bin produced by \code{\link{bin_prop}}.
#' This function is helpful for determining the obligatory parameter input to 
#' \code{\link{transform_to_elogit}}.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{bin_prop}}.
#' @return A printed summary of the number of samples in each bin.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Determine the number of samples per bin...
#' check_samples_per_bin(dat)
#' }
check_samples_per_bin <- function (data = data) {
  samples <- max(data$IA_0_C)
  rate <- abs(data$Time[2] - data$Time[1])
  print(paste("There are", samples, "samples in each bin."))
  print(paste("One data point every", rate, "millisecond(s)"))
}


#' Determine the sampling rate present in the data 
#' 
#' \code{check_samplingrate} determines the sampling rate in the data.
#' This function is helpful for determining the obligatory parameter input to 
#' \code{\link{bin_prop}}. If different sampling rates were used, the 
#' function adds a sampling rate column, which can be used to subset the 
#' data for further processing.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{select_recorded_eye}}.
#' @param ReturnData A logical indicating whether to return a data table containing 
#' a new column called SamplingRate
#' @return A printed summary and/or a data table object 
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Determine the sampling rate...
#' check_samplingrate(dat)
#' }
check_samplingrate <- function(data = data, ReturnData = FALSE) {
  ReturnData <- ReturnData
  
  tmp <- data %>%
    group_by(Event) %>%
    mutate(., SamplingRate = 1000 / (Time[2] - Time[1]))
  print(paste("Sampling rate(s) present in the data are:", unique(tmp$SamplingRate), "Hz."))
  
  if (length(unique(tmp$SamplingRate)) > 1) {
    warning("There are multiple sampling rates present in the data. Please use the ReturnData parameter to include a sampling rate column in the dataset. This can be used to subset by sampling rate before proceeding with the remaining preprocessing operations.")
  } 
  
  if (ReturnData == TRUE) {
  return(tmp)
  }
  
}



#' Determine downsampling options based on current sampling rate
#' 
#' \code{ds_options} determines the possible rates to which 
#' the current sampling rate can downsampled. It then prints the 
#' options in both bin size (milliseconds) and corresponding 
#' sampling rate (Hertz).
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' 
#' @param SamplingRate A positive integer indicating the sampling rate (in Hertz) 
#' used to record the gaze data, which can be determined with the function 
#' \code{\link{check_samplingrate}}.
#' @return A printed summary of options (bin size and rate) for downsampling.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Determine downsampling options...
#' ds_options(SamplingRate = 1000)
#' }
ds_options <- function(SamplingRate=SamplingRate) {
  SamplingRate = SamplingRate
  for (x in 1:100) {
  if (x %% (1000/SamplingRate) == 0) {
    if ((1000/x) %% 1 == 0) {
    print(paste("Bin size:", x, "ms;", "Downsampled rate:", 1000/x, "Hz"))
    }
  }
  }
}




#' Check which eyes were recorded during the experiment
#' 
#' \code{check_eye_recording} quickly checks if the dataset contains gaze data
#' in both the Right and Left interest area columns. It prints a summary and 
#' suggests which setting to use for the \code{Recording} parameter in the 
#' function \code{\link{select_recorded_eye}}.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{create_time_series}}.
#' @return Text feedback and instruction.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Create a unified columns for the gaze data...
#' check_eye_recording(dat)
#' }
check_eye_recording <- function(data = data) {
  
  if (sum(data$LEFT_INTEREST_AREA_ID) > 0 & sum(data$RIGHT_INTEREST_AREA_ID) > 0) {
    print("The dataset contains recordings for both eyes. If any participants had both eyes tracked, set the Recording parameter in select_recorded_eye() to 'LandR'. If participants had either the left OR the right eye tracked, set the Recording parameter in select_recorded_eye() to 'LorR'.")
  } else if (sum(data$LEFT_INTEREST_AREA_ID) > 0 & sum(data$RIGHT_INTEREST_AREA_ID) == 0) {
    print("The dataset contains recordings for ONLY the left eye. Set the Recording parameter in select_recorded_eye() to 'L'.")
  } else if (sum(data$LEFT_INTEREST_AREA_ID) == 0 & sum(data$RIGHT_INTEREST_AREA_ID) > 0) {
    print("The dataset contains recordings for ONLY the right eye. Set the Recording parameter in select_recorded_eye() to 'R'.")
  }
  
}




#' Rename default column names for interest areas.
#' 
#' \code{rename_columns} will replace the default numerical coding of the 
#' interest area columns with more meaningful user-specified names. For example,
#' IA_1_C and IA_1_P could be converted to IA_Target_C and IA_Target_P. Again, 
#' this will work for upto 8 interest areas.
#' 
#' @export
#' @import dplyr
#' @import tidyr
#' @import lazyeval
#' 
#' @param data A data table object output by either \code{\link{bin_prop}}. 
#' \code{\link{transform_to_elogit}}, or \code{\link{create_binomial}}.
#' @param Labels A named character vector specifying the interest areas and the
#' desired names to be inserted in place of the numerical labelling.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # For renaming default interest area columns
#' dat2 <- rename_columns(dat, Labels = c(IA1="Target", IA2="Rhyme", 
#'                            IA3="OnsetComp", IA4="Distractor")) 
#' }
rename_columns <- function(data = data, Labels = Labels) {
  
  Labels <- Labels
  tmp <- data

  if (length(names(Labels))>8) {
    stop("You have more than 8 interest areas.")
  } else {
    print(paste("Renaming", length(names(Labels)), "interest areas.", sep = " "))
  }
    
  Labels <- c("0" = "outside", Labels)
  
  NoIA <- length(names(Labels))
  
  for (x in 1:NoIA) {
    Labels[[x]] <- paste("_",Labels[[x]],"_", sep = "")
    names(Labels)[x] <- paste("_",x-1,"_", sep = "")
    tmp<-setNames(tmp, gsub(names(Labels)[x],Labels[[x]],names(tmp)))
  }
  
  return(tmp)
  
}

