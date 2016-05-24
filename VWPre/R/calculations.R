#' Bins the sample data and calculates proportion looks by interest area
#' 
#' \code{bin_prop} calculates the proportion of looks (samples) to each 
#' interest area in a particular window of time (bin size). This function first
#' checks to see if the procedure is possible given the sampling rate and 
#' desired bin size. It then performs the calculation and downsampling, returning
#' new columns corresponding to each interest area ID (e.g., 'IA_1_C', 'IA_1_P').
#' The extention '_c' indicates the count of samples in the bin and the 
#' extension '_P' indicates the proportion. N.B.: This function will work for 
#' data with a maximum of 8 interest areas.
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{select_recorded_eye}}
#' or \code{\link{check_samplingrate}}.
#' @param NoIA A positive integer indicating the number of interest areas defined 
#' when creating the study. 
#' @param BinSize A positive integer indicating the size of the binning window 
#' (in milliseconds).
#' @param SamplingRate A positive integer indicating the sampling rate (in Hertz) 
#' used to record the gaze data, which can be determined with the function 
#' \code{\link{check_samplingrate}}.
#' @return A data table with additional columns (the number of which depends on 
#' the number of interest areas specified) added to \code{data}.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Bin samples and calculation proportions...
#' df <- bin_prop(dat, NoIA = 4, BinSize = 20, SamplingRate = 1000)
#' }
bin_prop <- function(data = data, NoIA = NoIA, BinSize = BinSize, SamplingRate = SamplingRate) {
  NoIA = NoIA
  samplerate <-  data$Time[2] - data$Time[1]
  BinSize = BinSize
  SamplesPerBin <- (SamplingRate / 1000) * BinSize
  if (BinSize %% samplerate != 0) {
    stop("Sample frequency of data is not a multiple of the target frequency.")
  } else {
    print("Sampling rate OK. You're good to go!")
  }
  data$DS <- data$Time %/% BinSize
  data$DS <- data$DS * BinSize 
  if (NoIA == 0) {
    stop("You must have at least one interest area!")
  } else if (NoIA > 8) {
    stop("You have more than 8 interest areas; you must modify this function.")
  } else if (NoIA == 1) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin)
  } else if (NoIA == 2) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin)
  } else if (NoIA == 3) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin)
  } else if (NoIA == 4) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)]),
             IA_4_C = length(IA_ID[which(IA_ID == 4)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin,
             IA_4_P = IA_4_C / SamplesPerBin)
  } else if (NoIA == 5) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)]),
             IA_4_C = length(IA_ID[which(IA_ID == 4)]),
             IA_5_C = length(IA_ID[which(IA_ID == 5)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin,
             IA_4_P = IA_4_C / SamplesPerBin,
             IA_5_P = IA_5_C / SamplesPerBin)
  } else if (NoIA == 6) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)]),
             IA_4_C = length(IA_ID[which(IA_ID == 4)]),
             IA_5_C = length(IA_ID[which(IA_ID == 5)]),
             IA_6_C = length(IA_ID[which(IA_ID == 6)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin,
             IA_4_P = IA_4_C / SamplesPerBin,
             IA_5_P = IA_5_C / SamplesPerBin,
             IA_6_P = IA_6_C / SamplesPerBin)
  } else if (NoIA == 7) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)]),
             IA_4_C = length(IA_ID[which(IA_ID == 4)]),
             IA_5_C = length(IA_ID[which(IA_ID == 5)]),
             IA_6_C = length(IA_ID[which(IA_ID == 6)]),
             IA_7_C = length(IA_ID[which(IA_ID == 7)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin,
             IA_4_P = IA_4_C / SamplesPerBin,
             IA_5_P = IA_5_C / SamplesPerBin,
             IA_6_P = IA_6_C / SamplesPerBin,
             IA_7_P = IA_7_C / SamplesPerBin)
  } else if (NoIA == 8) {
    data %>%
      group_by(Event, DS) %>%
      mutate(., IA_0_C = length(IA_ID[which(IA_ID == 0)]),
             IA_1_C = length(IA_ID[which(IA_ID == 1)]),
             IA_2_C = length(IA_ID[which(IA_ID == 2)]),
             IA_3_C = length(IA_ID[which(IA_ID == 3)]),
             IA_4_C = length(IA_ID[which(IA_ID == 4)]),
             IA_5_C = length(IA_ID[which(IA_ID == 5)]),
             IA_6_C = length(IA_ID[which(IA_ID == 6)]),
             IA_7_C = length(IA_ID[which(IA_ID == 7)]),
             IA_8_C = length(IA_ID[which(IA_ID == 8)])) %>% 
      filter(., Time %in% DS) %>%
      mutate(., IA_0_P = IA_0_C / SamplesPerBin,
             IA_1_P = IA_1_C / SamplesPerBin,
             IA_2_P = IA_2_C / SamplesPerBin,
             IA_3_P = IA_3_C / SamplesPerBin,
             IA_4_P = IA_4_C / SamplesPerBin,
             IA_5_P = IA_5_C / SamplesPerBin,
             IA_6_P = IA_6_C / SamplesPerBin,
             IA_7_P = IA_7_C / SamplesPerBin,
             IA_8_P = IA_8_C / SamplesPerBin)
  }
  
}




#' Transforms proportion looks to empirical logits.
#' 
#' \code{transform_to_elogit} transforms the proportion of looks for 
#' each interest area to empirical logits. Proportions are inherently bound 
#' between 0 and 1 and are therefore not suitable for some types of analysis. 
#' Logits provide an unbounded measure, though range from negative infinity to 
#' infinity, so it is important to know that this logit function adds a constant 
#' (hence, empirical logit). Additionally this calculates weights which estimate 
#' the variance in each bin (because the variance of the logit depends on the 
#' mean). This is important for regression analyses. N.B.: This function will 
#' work for data with a maximum of 8 interest areas.
#' 
#' These calculations are taken from:
#' Barr, D. J., (2008) Analyzing 'visual world' eyetracking data using 
#' multilevel logistic regression, \emph{Journal of Memory and Language}, 
#' \emph{59}(4), 457--474.
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' 
#' @param data A data table object output by \code{\link{select_recorded_eye}}.
#' @param NoIA A positive integer indicating the number of interest areas defined 
#' when creating the study. 
#' @param SamplesPerBin A positive integer indicating the number of samples in
#' each bin, which can be determines with \code{\link{check_samples_per_bin}}.
#' @param Constant A positive number used for the empirical logit and weights
#' calculation; by default, 0.5 as in Barr (2008).
#' @return A data table with additional columns (the number of which depends on 
#' the number of interest areas specified) added to \code{data}.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Convert proportions to empirical logits and calculate weights...
#' df <- transform_to_elogit(dat, NoIA = 4, SamplesPerBin = 20, Constant = 0.5)
#' }
transform_to_elogit <- function(data = data, NoIA = NoIA, SamplesPerBin = SamplesPerBin,
                                Constant = 0.5) {
  NoIA = NoIA
  SamplesPerBin = SamplesPerBin
  Constant = Constant
  
  if (NoIA > 8) {
    stop("You have too many interest areas. You must modify the function.")
  }
  else if (NoIA == 1) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant))
  }
  else if (NoIA == 2) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant))
  }
  else if (NoIA == 3) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant))
  }
  else if (NoIA == 4) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_4_ELogit = log((IA_4_C + Constant) / (SamplesPerBin - IA_4_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant),
             IA_4_wts = 1/(IA_4_C + Constant) + 1/(SamplesPerBin - IA_4_C + Constant))
  }
  else if (NoIA == 5) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_4_ELogit = log((IA_4_C + Constant) / (SamplesPerBin - IA_4_C + Constant)),
             IA_5_ELogit = log((IA_5_C + Constant) / (SamplesPerBin - IA_5_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant),
             IA_4_wts = 1/(IA_4_C + Constant) + 1/(SamplesPerBin - IA_4_C + Constant),
             IA_5_wts = 1/(IA_5_C + Constant) + 1/(SamplesPerBin - IA_5_C + Constant))
  }
  else if (NoIA == 6) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_4_ELogit = log((IA_4_C + Constant) / (SamplesPerBin - IA_4_C + Constant)),
             IA_5_ELogit = log((IA_5_C + Constant) / (SamplesPerBin - IA_5_C + Constant)),
             IA_6_ELogit = log((IA_6_C + Constant) / (SamplesPerBin - IA_6_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant),
             IA_4_wts = 1/(IA_4_C + Constant) + 1/(SamplesPerBin - IA_4_C + Constant),
             IA_5_wts = 1/(IA_5_C + Constant) + 1/(SamplesPerBin - IA_5_C + Constant),
             IA_6_wts = 1/(IA_6_C + Constant) + 1/(SamplesPerBin - IA_6_C + Constant))
  }
  else if (NoIA == 7) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_4_ELogit = log((IA_4_C + Constant) / (SamplesPerBin - IA_4_C + Constant)),
             IA_5_ELogit = log((IA_5_C + Constant) / (SamplesPerBin - IA_5_C + Constant)),
             IA_6_ELogit = log((IA_6_C + Constant) / (SamplesPerBin - IA_6_C + Constant)),
             IA_7_ELogit = log((IA_7_C + Constant) / (SamplesPerBin - IA_7_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant),
             IA_4_wts = 1/(IA_4_C + Constant) + 1/(SamplesPerBin - IA_4_C + Constant),
             IA_5_wts = 1/(IA_5_C + Constant) + 1/(SamplesPerBin - IA_5_C + Constant),
             IA_6_wts = 1/(IA_6_C + Constant) + 1/(SamplesPerBin - IA_6_C + Constant),
             IA_7_wts = 1/(IA_7_C + Constant) + 1/(SamplesPerBin - IA_7_C + Constant))
  }
  else if (NoIA == 8) {
    data %>%
      mutate(., IA_0_ELogit = log((IA_0_C + Constant) / (SamplesPerBin - IA_0_C + Constant)),
             IA_1_ELogit = log((IA_1_C + Constant) / (SamplesPerBin - IA_1_C + Constant)),
             IA_2_ELogit = log((IA_2_C + Constant) / (SamplesPerBin - IA_2_C + Constant)),
             IA_3_ELogit = log((IA_3_C + Constant) / (SamplesPerBin - IA_3_C + Constant)),
             IA_4_ELogit = log((IA_4_C + Constant) / (SamplesPerBin - IA_4_C + Constant)),
             IA_5_ELogit = log((IA_5_C + Constant) / (SamplesPerBin - IA_5_C + Constant)),
             IA_6_ELogit = log((IA_6_C + Constant) / (SamplesPerBin - IA_6_C + Constant)),
             IA_7_ELogit = log((IA_7_C + Constant) / (SamplesPerBin - IA_7_C + Constant)),
             IA_8_ELogit = log((IA_8_C + Constant) / (SamplesPerBin - IA_8_C + Constant)),
             IA_0_wts = 1/(IA_0_C + Constant) + 1/(SamplesPerBin - IA_0_C + Constant),
             IA_1_wts = 1/(IA_1_C + Constant) + 1/(SamplesPerBin - IA_1_C + Constant),
             IA_2_wts = 1/(IA_2_C + Constant) + 1/(SamplesPerBin - IA_2_C + Constant),
             IA_3_wts = 1/(IA_3_C + Constant) + 1/(SamplesPerBin - IA_3_C + Constant),
             IA_4_wts = 1/(IA_4_C + Constant) + 1/(SamplesPerBin - IA_4_C + Constant),
             IA_5_wts = 1/(IA_5_C + Constant) + 1/(SamplesPerBin - IA_5_C + Constant),
             IA_6_wts = 1/(IA_6_C + Constant) + 1/(SamplesPerBin - IA_6_C + Constant),
             IA_7_wts = 1/(IA_7_C + Constant) + 1/(SamplesPerBin - IA_7_C + Constant),
             IA_8_wts = 1/(IA_8_C + Constant) + 1/(SamplesPerBin - IA_8_C + Constant))
  }
}





#' Creates a success/failure column for each IA based on counts.
#' 
#' \code{create_binomial} uses interest area count columns to create 
#' a success/failure column for each IA which is suitable for logistic regression. 
#' N.B.: This function will work for data with a maximum of 8 interest areas.
#' 
#' @export
#' @import dplyr
#' @import lazyeval
#' 
#' @param data A data table object output by either \code{\link{bin_prop}} or
#' \code{\link{transform_to_elogit}}.
#' @param NoIA A positive integer indicating the number of interest areas defined 
#' when creating the study. 
#' @return A data table with additional columns (the number of which depends on 
#' the number of interest areas specified) added to \code{data}.
#' @examples
#' \dontrun{
#' library(VWPre)
#' # Create binomial (success/failure) column...
#' df <- create_binomial(data = dat, NoIA = 4)
#' }
create_binomial <- function(data = data, NoIA = NoIA) {
  
  data <- data %>% ungroup()
  
  if (NoIA > 8) {
    stop("You have too many interest areas. You must modify the function.")
  }
  else if (NoIA == 1) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C,
               IA_1_off = IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off)
    return(tmp)
  }
  
  else if (NoIA == 2) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C,
               IA_1_off = IA_0_C + IA_2_C ,
               IA_2_off = IA_1_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off)
    return(tmp)
  }
  
  else if (NoIA == 3) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C + IA_3_C,
               IA_1_off = IA_0_C + IA_2_C + IA_3_C,
               IA_2_off = IA_1_C + IA_0_C + IA_3_C,
               IA_3_off = IA_1_C + IA_2_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off)
    return(tmp)
  }
  
  else if (NoIA == 4) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C,
               IA_1_off = IA_0_C + IA_2_C + IA_3_C + IA_4_C,
               IA_2_off = IA_1_C + IA_0_C + IA_3_C + IA_4_C,
               IA_3_off = IA_1_C + IA_2_C + IA_0_C + IA_4_C,
               IA_4_off = IA_1_C + IA_2_C + IA_3_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    tmp$IA_4_Looks = cbind(tmp$IA_4_C, tmp$IA_4_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off, -IA_4_off)
    return(tmp)
  }
  
  else if (NoIA == 5) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C,
               IA_1_off = IA_0_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C,
               IA_2_off = IA_1_C + IA_0_C + IA_3_C + IA_4_C + IA_5_C,
               IA_3_off = IA_1_C + IA_2_C + IA_0_C + IA_4_C + IA_5_C,
               IA_4_off = IA_1_C + IA_2_C + IA_3_C + IA_0_C + IA_5_C,
               IA_5_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    tmp$IA_4_Looks = cbind(tmp$IA_4_C, tmp$IA_4_off)
    tmp$IA_5_Looks = cbind(tmp$IA_5_C, tmp$IA_5_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off, -IA_4_off, -IA_5_off)
    return(tmp)
  }
  
  else if (NoIA == 6) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C,
               IA_1_off = IA_0_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C,
               IA_2_off = IA_1_C + IA_0_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C,
               IA_3_off = IA_1_C + IA_2_C + IA_0_C + IA_4_C + IA_5_C + IA_6_C,
               IA_4_off = IA_1_C + IA_2_C + IA_3_C + IA_0_C + IA_5_C + IA_6_C,
               IA_5_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_0_C + IA_6_C,
               IA_6_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    tmp$IA_4_Looks = cbind(tmp$IA_4_C, tmp$IA_4_off)
    tmp$IA_5_Looks = cbind(tmp$IA_5_C, tmp$IA_5_off)
    tmp$IA_6_Looks = cbind(tmp$IA_6_C, tmp$IA_6_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off, -IA_4_off, -IA_5_off, -IA_6_off)
    return(tmp)
  }
  
  else if (NoIA == 7) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C + IA_7_C,
               IA_1_off = IA_0_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C + IA_7_C,
               IA_2_off = IA_1_C + IA_0_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C + IA_7_C,
               IA_3_off = IA_1_C + IA_2_C + IA_0_C + IA_4_C + IA_5_C + IA_6_C + IA_7_C,
               IA_4_off = IA_1_C + IA_2_C + IA_3_C + IA_0_C + IA_5_C + IA_6_C + IA_7_C,
               IA_5_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_0_C + IA_6_C + IA_7_C,
               IA_6_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_0_C + IA_7_C,
               IA_7_off = IA_1_C + IA_2_C + IA_3_C + IA_4_C + IA_5_C + IA_6_C + IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    tmp$IA_4_Looks = cbind(tmp$IA_4_C, tmp$IA_4_off)
    tmp$IA_5_Looks = cbind(tmp$IA_5_C, tmp$IA_5_off)
    tmp$IA_6_Looks = cbind(tmp$IA_6_C, tmp$IA_6_off)
    tmp$IA_7_Looks = cbind(tmp$IA_7_C, tmp$IA_7_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off, -IA_4_off, -IA_5_off, -IA_6_off, -IA_7_off)
    return(tmp)
  }
  
  else if (NoIA == 8) {
    tmp <- data %>% 
      group_by(Event) %>%
      do(
        mutate(., IA_0_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_1_off = IA_0_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_2_off = IA_1_C  +  IA_0_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_3_off = IA_1_C  +  IA_2_C  +  IA_0_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_4_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_0_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_5_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_0_C  +  IA_6_C  +  IA_7_C  +  IA_8_C,
               IA_6_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_0_C  +  IA_7_C  +  IA_8_C,
               IA_7_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_0_C  +  IA_8_C,
               IA_8_off = IA_1_C  +  IA_2_C  +  IA_3_C  +  IA_4_C  +  IA_5_C  +  IA_6_C  +  IA_7_C  +  IA_0_C)
      ) 
    
    tmp$IA_0_Looks = cbind(tmp$IA_0_C, tmp$IA_0_off)
    tmp$IA_1_Looks = cbind(tmp$IA_1_C, tmp$IA_1_off)
    tmp$IA_2_Looks = cbind(tmp$IA_2_C, tmp$IA_2_off)
    tmp$IA_3_Looks = cbind(tmp$IA_3_C, tmp$IA_3_off)
    tmp$IA_4_Looks = cbind(tmp$IA_4_C, tmp$IA_4_off)
    tmp$IA_5_Looks = cbind(tmp$IA_5_C, tmp$IA_5_off)
    tmp$IA_6_Looks = cbind(tmp$IA_6_C, tmp$IA_6_off)
    tmp$IA_7_Looks = cbind(tmp$IA_7_C, tmp$IA_7_off)
    tmp$IA_8_Looks = cbind(tmp$IA_8_C, tmp$IA_8_off)
    
    tmp <- select(tmp, -IA_0_off, -IA_1_off, -IA_2_off, -IA_3_off, -IA_4_off, -IA_5_off, -IA_6_off, -IA_7_off, -IA_8_off)
    return(tmp)
  }
  
}

