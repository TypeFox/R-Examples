


#' Biascorrect the input timeseries or hyfo dataset
#' 
#' Biascorrect the input time series or dataset, the input time series or dataset should consist of observation, hindcast, and forecast.
#' observation and hindcast should belong to the same period, in order to calibrate. Then the modified forecast
#' will be returned. If the input is a time series, first column should be date column and rest columns should be 
#' the value column. If the input is a hyfo dataset, the dataset should be the result of \code{loadNcdf}, or a list
#' file with the same format. 
#' 
#' 
#' @param frc a hyfo grid data output or a dataframe (time series) consists of Date column and one or more value columns, 
#' representing the forecast to be calibrated.
#' @param hindcast a hyfo grid data output or a dataframe(time series) consists of Date column and one or more value columns, 
#' representing the hindcast data. This data will be used in the calibration of the forecast, so it's better to have the same date period as
#' observation data. Check details for more information.
#' @param obs a hyfo grid data output or a dataframe (time series) consists of Date column and one or more value columns, 
#' representing the observation data.
#' @param method bias correct method, including 'delta', 'scaling'..., default is 'scaling'
#' @param scaleType only when the method "scaling" is chosen, scaleType will be available. Two different types
#' of scaling method, 'add' and 'multi', which means additive and multiplicative scaling method. More info check 
#' details. Default scaleType is 'multi'.
#' @param preci If the precipitation is biascorrected, then you have to assign \code{preci = TRUE}. Since for
#' precipitation, some biascorrect methods may not apply to, or some methods are specially for precipitation. 
#' Default is FALSE, refer to details.
#' @param prThreshold The minimum value that is considered as a non-zero precipitation. Default to 1 (assuming mm).
#' @param extrapolate When use 'eqm' method, and 'no' is set, modified frc is bounded by the range of obs.
#' If 'constant' is set, modified frc is not bounded by the range of obs. Default is 'no'.
#' @details 
#' 
#' Since climate forecast is based on global condition, when downscaling to different regions, it may include
#' some bias, biascorrection is used then to fix the bias.
#' 
#' \strong{Hindcast}
#' 
#' In order to bias correct, we need to pick up some data from the forecast to train with
#' the observation, which is called hindcast in this function. Using hindcast and observation, 
#' the program can analyze the bias and correct the bias in the forecast. 
#' 
#' Hindcast should have \strong{EVERY} attributes that forecast has.
#' 
#' Hindcast is also called re-forecast, is the forecast of the past. E.g. you have a forecast from year 2000-2010, assuming now you are in 2005. So from 2000-2005, this period
#' is the hindcast period, and 2005-2010, this period is the forecast period.
#'
#' Hindcast can be the same as forecast, i.e., you can use forecast itself as hindcast to train the bias correction.
#'
#'
#' \strong{How it works}
#' 
#' Forecast product has to be calibrated, usually the system is doing forecast in real time. So, e.g., if the 
#' forecast starts from year 2000, assuming you are in year 2003, then you will have 3 years' hindcast 
#' data (year 2000-2003), which can be used to calibrate. And your forecast period is (2003-2004)
#' 
#' E.g. you have observation from 2001-2002, this is your input obs. Then you can take the same 
#' period (2001-2002) from the forecast, which is the hindcast period. For forecast, you can take any period.
#' The program will evaluate the obs and hindcast, to get the modification of the forecast, and then add the 
#' modification to the forecast data.
#' 
#' The more categorized input, the more accurate result you will get. E.g., if you want to 
#' bias correct a forecast for winter season. So you'd better to extract all the winter period
#' in the hindcast and observation to train. \code{extractPeriod} can be used for this purpose.
#' 
#' \strong{method}
#' 
#' Different methods used in the bias correction. Among which, delta, scaling can be applied
#' to different kinds of parameters, with no need to set \code{preci}; eqm has two conditions for rainfall data and other data,
#' it needs user to input \code{preci = TRUE/FALSE} to point to different conditions; gqm is
#' designed for rainfall data, so \code{preci = TRUE} needs to be set.
#' 
#' \strong{delta}
#' 
#' This method consists on adding to the observations the mean change signal (delta method). 
#' This method is applicable to any kind of variable but it is preferable to avoid it for bounded variables
#'  (e.g. precipitation, wind speed, etc.) because values out of the variable range could be obtained 
#'  (e.g. negative wind speeds...)
#'  
#'  \strong{scaling}
#'  
#' This method consists on scaling the simulation  with the difference (additive) or quotient (multiplicative) 
#' between the observed and simulated means in the train period. The \code{additive} or \code{multiplicative}
#' correction is defined by parameter \code{scaling.type} (default is \code{additive}).
#' The additive version is preferably applicable to unbounded variables (e.g. temperature) 
#' and the multiplicative to variables with a lower bound (e.g. precipitation, because it also preserves the frequency). 
#'  
#'  \strong{eqm}
#'  
#' Empirical Quantile Mapping. This is a very extended bias correction method which consists on calibrating the simulated Cumulative Distribution Function (CDF) 
#' by adding to the observed quantiles both the mean delta change and the individual delta changes in the corresponding quantiles. 
#' This method is applicable to any kind of variable.
#' 
#' It can keep the extreme value, if you choose constant extrapolation method. But then you will face the risk
#' that the extreme value is an error.
#'  
#'  \strong{gqm}
#'  
#' Gamma Quantile Mapping. This method is described in Piani et al. 2010 and is applicable only to precipitation. It is based on the initial assumption that both observed
#' and simulated intensity distributions are well approximated by the gamma distribution, therefore is a parametric q-q map 
#' that uses the theorical instead of the empirical distribution. 
#'  
#' It can somehow filter some extreme values caused by errors, while keep the extreme value. Seems more reasonable.
#' Better have a long period of training, and the if the forecast system is relatively stable.
#' 
#' It is a generic function, if in your case you need to debug, please see \code{?debug()} 
#' for how to debug S4 method.
#'  
#' @examples 
#' 
#' ######## hyfo grid file biascorrection
#' ########
#' 
#' # If your input is obtained by \code{loadNcdf}, you can also directly biascorrect
#' # the file.
#' 
#' # First load ncdf file.
#' filePath <- system.file("extdata", "tnc.nc", package = "hyfo")
#' varname <- getNcdfVar(filePath)    
#' nc <- loadNcdf(filePath, varname)
#' 
#' data(tgridData)
#' # Since the example data, has some NA values, the process will include some warning #message, 
#' # which can be ignored in this case.
#' 
#' 
#' 
#' 
#' # Then we will use nc data as forecasting data, and use itself as hindcast data,
#' # use tgridData as observation.
#' newFrc <- biasCorrect(nc, nc, tgridData)  
#' newFrc <- biasCorrect(nc, nc, tgridData, scaleType = 'add')   
#' newFrc <- biasCorrect(nc, nc, tgridData, method = 'eqm', extrapolate = 'constant', 
#' preci = TRUE) 
#' newFrc <- biasCorrect(nc, nc, tgridData, method = 'gqm', preci = TRUE) 
#' 
#' 
#' ######## Time series biascorrection
#' ########
#' 
#' # Use the time series from testdl as an example, we take frc, hindcast and obs from testdl.
#' data(testdl)
#' 
#' # common period has to be extracted in order to better train the forecast.
#' 
#' datalist <- extractPeriod(testdl, startDate = '1994-1-1', endDate = '1995-10-1')
#' 
#' frc <- datalist[[1]]
#' hindcast <- datalist[[2]]
#' obs <- datalist[[3]]
#' 
#' 
#' # The data used here is just for example, so there could be negative data.
#' 
#' # default method is scaling, with 'multi' scaleType
#' frc_new <- biasCorrect(frc, hindcast, obs)
#' 
#' # for precipitation data, extra process needs to be executed, so you have to tell
#' # the program that it is a precipitation data.
#' 
#' frc_new1 <- biasCorrect(frc, hindcast, obs, preci = TRUE)
#' 
#' # You can use other scaling methods to biascorrect.
#' frc_new2 <- biasCorrect(frc, hindcast, obs, scaleType = 'add')
#' 
#' # 
#' frc_new3 <- biasCorrect(frc, hindcast, obs, method = 'eqm', preci = TRUE)
#' frc_new4 <- biasCorrect(frc, hindcast, obs, method = 'gqm', preci = TRUE)
#' 
#' plotTS(obs, frc, frc_new, frc_new1, frc_new2, frc_new3, frc_new4, plot = 'cum')
#' 
#' # You can also give name to this input list.
#' TSlist <- list(obs, frc, frc_new, frc_new1, frc_new2, frc_new3, frc_new4)
#' names(TSlist) <- c('obs', 'frc', 'delta', 'delta_preci', 'scale', 'eqm', 'gqm')
#' plotTS(list = TSlist, plot = 'cum')
#' 
#' 
#' # If the forecasts you extracted only has incontinuous data for certain months and years, e.g.,
#' # for seasonal forecasting, forecasts only provide 3-6 months data, so the case can be 
#' # for example Dec, Jan and Feb of every year from year 1999-2005.
#' # In such case, you need to extract certain months and years from observed time series.
#' # extractPeriod() can be then used.
#'   
#'   
#'
#'
#'
#' # More examples can be found in the user manual on http://yuanchao-xu.github.io/hyfo/
#' 
#' 
#' @references 
#' Bias correction methods come from \code{biasCorrection} from \code{dowscaleR}
#' 
#' \itemize{
#' 
#' \item Santander Meteorology Group (2015). downscaleR: Climate data manipulation and statistical downscaling. R
#' package version 0.6-0. https://github.com/SantanderMetGroup/downscaleR/wiki
#' 
#' \item R.A.I. Wilcke, T. Mendlik and A. Gobiet (2013) Multi-variable error correction of regional climate models. Climatic Change, 120, 871-887
#' 
#' \item A. Amengual, V. Homar, R. Romero, S. Alonso, and C. Ramis (2012) A Statistical Adjustment of Regional Climate Model Outputs to Local Scales: Application to Platja de Palma, Spain. J. Clim., 25, 939-957
#' 
#' \item C. Piani, J. O. Haerter and E. Coppola (2009) Statistical bias correction for daily precipitation in regional climate models over Europe, Theoretical and Applied Climatology, 99, 187-192
#' 
#' \item O. Gutjahr and G. Heinemann (2013) Comparing precipitation bias correction methods for high-resolution regional climate simulations using COSMO-CLM, Theoretical and Applied Climatology, 114, 511-529
#' }
#' 
#' @author Yuanchao Xu \email{xuyuanchao37@@gmail.com }
#' @importFrom methods setMethod
#' @export
#' 
setGeneric('biasCorrect', function(frc, hindcast, obs, method = 'scaling', scaleType = 'multi', 
                                   preci = FALSE, prThreshold = 0, extrapolate = 'no') {
  standardGeneric('biasCorrect')
})

#' @describeIn biasCorrect
setMethod('biasCorrect', signature('data.frame', 'data.frame', 'data.frame'),
          function(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
            result <- biasCorrect.TS(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate)
            return(result)
          })

#' @describeIn biasCorrect
setMethod('biasCorrect', signature('list', 'list', 'list'), 
          function(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
            result <- biasCorrect.list(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate)
            return(result)
          })


biasCorrect.TS <- function(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
  # First check if the first column is Date
  if (!grepl('-|/', obs[1, 1]) | !grepl('-|/', hindcast[1, 1]) | !grepl('-|/', frc[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
         and use as.Date to convert.If your input is a hyfo dataset, put input = "hyfo" as an
         argument, check help for more info.')
  }
  # change to date type is easier, but in case in future the flood part is added, Date type doesn't have
  # hour, min and sec, so, it's better to convert it into POSIxlt.
  
  # if condition only accepts one condition, for list comparison, there are a lot of conditions, better
  # further process it, like using any.
  if (any(as.POSIXlt(hindcast[, 1]) != as.POSIXlt(obs[, 1]))) {
    warning('time of obs and time of hindcast are not the same, which may cause inaccuracy in 
                      the calibration.')
  }
  n <- ncol(frc)
  
  # For every column, it's biascorrected respectively.
  frc_data <- lapply(2:n, function(x) biasCorrect_core(frc[, x], hindcast[, x], obs[, 2], method = method,
                                                       scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                       extrapolate = extrapolate))
  frc_data <- do.call('cbind', frc_data)
  rownames(frc_data) <- NULL
  
  names <- colnames(frc)
  frc_new <- data.frame(frc[, 1], frc_data)
  colnames(frc_new) <- names
  return(frc_new)
}

biasCorrect.list <- function(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
  ## Check if the data is a hyfo grid data.
  checkHyfo(frc, hindcast, obs)
  
  hindcastData <- hindcast$Data
  obsData <- obs$Data
  frcData <- frc$Data
  
  ## save frc dimension order, at last, set the dimension back to original dimension
  frcDim <- attributes(frcData)$dimensions
  
  ## ajust the dimension into general dimension order.
  hindcastData <- adjustDim(hindcastData, ref = c('lon', 'lat', 'time'))
  obsData <- adjustDim(obsData, ref = c('lon', 'lat', 'time'))
  
  ## CheckDimLength, check if all the input dataset has different dimension length
  # i.e. if they all have the same lon and lat number.
  checkDimLength(frcData, hindcastData, obsData, dim = c('lon', 'lat'))
  
  
  # Now real bias correction is executed.
  
  memberIndex <- grepAndMatch('member', attributes(frcData)$dimensions)
  
  # For dataset that has a member part 
  if (length(memberIndex) != 0) {
    # check if frcData and hindcastData has the same dimension and length.
    checkDimLength(frcData, hindcastData, dim = 'member')
    
    frcData <- adjustDim(frcData, ref = c('lon', 'lat', 'time', 'member'))
    
    # The following code may speed up because it doesn't use for loop.
    # It firstly combine different array into one array. combine the time 
    # dimension of frc, hindcast and obs. Then use apply, each time extract 
    # the total time dimension, and first part is frc, second is hindcast, third
    # is obs. Then use these three parts to bias correct. All above can be written
    # in one function and called within apply. But too complicated to understand,
    # So save it for future use maybe.
    
    #       for (member in 1:dim(frcData)[4]) {
    #         totalArr <- array(c(frcData[,,, member], hindcastData[,,, member], obsData),
    #                           dim = c(dim(frcData)[1], dim(frcData)[2], 
    #                                                        dim(frcData)[3] + dim(hindcastData)[3] + dim(obsData)[3]))
    #       }
    
    
    for (member in 1:dim(frcData)[4]) {
      for (lon in 1:dim(frcData)[1]) {
        for (lat in 1:dim(frcData)[2]) {
          frcData[lon, lat,, member] <- biasCorrect_core(frcData[lon, lat,,member], hindcastData[lon, lat,, member], obsData[lon, lat,], method = method,
                                                         scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                         extrapolate = extrapolate)
        }
      }
    }
  } else {
    frcData <- adjustDim(frcData, ref = c('lon', 'lat', 'time'))
    for (lon in 1:dim(frcData)[1]) {
      for (lat in 1:dim(frcData)[2]) {
        frcData[lon, lat,] <- biasCorrect_core(frcData[lon, lat,], hindcastData[lon, lat,], obsData[lon, lat,], method = method,
                                               scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                               extrapolate = extrapolate)
      }
    }
  }
  
  frcData <- adjustDim(frcData, ref = frcDim)
  frc$Data <- frcData
  frc$biasCorrected_by <- method
  frc_new <- frc
  return(frc_new)
}






#' @importFrom MASS fitdistr
#' @importFrom stats ecdf quantile pgamma qgamma rgamma
#' 
#' @references 
#' Bias correction methods come from \code{biasCorrection} from \code{dowscaleR}
#' 
#' \itemize{
#' 
#' \item Santander Meteorology Group (2015). downscaleR: Climate data manipulation and statistical downscaling. R
#' package version 0.6-0. https://github.com/SantanderMetGroup/downscaleR/wiki
#' 
#' \item R.A.I. Wilcke, T. Mendlik and A. Gobiet (2013) Multi-variable error correction of regional climate models. Climatic Change, 120, 871-887
#' 
#' \item A. Amengual, V. Homar, R. Romero, S. Alonso, and C. Ramis (2012) A Statistical Adjustment of Regional Climate Model Outputs to Local Scales: Application to Platja de Palma, Spain. J. Clim., 25, 939-957
#' 
#' \item C. Piani, J. O. Haerter and E. Coppola (2009) Statistical bias correction for daily precipitation in regional climate models over Europe, Theoretical and Applied Climatology, 99, 187-192
#' 
#' \item O. Gutjahr and G. Heinemann (2013) Comparing precipitation bias correction methods for high-resolution regional climate simulations using COSMO-CLM, Theoretical and Applied Climatology, 114, 511-529
#' }
#' 
#' 
#' 
# this is only used to calculate the value column, 
biasCorrect_core <- function(frc, hindcast, obs, method, scaleType, preci, prThreshold, extrapolate){
  # If the variable is precipitation, some further process needs to be added.
  # The process is taken from downscaleR, to provide a more reasonable hindcast, used in the calibration.
  
  
  # check if frc, hindcast or obs are all na values
  if (!any(!is.na(obs)) | !any(!is.na(frc)) | !any(!is.na(hindcast))) {
    warning('In this cell, frc, hindcast or obs data is missing. No biasCorrection for this cell.')
    return(NA)
  }
  
  
  if (preci == TRUE) {
    preprocessHindcast_res <- preprocessHindcast(hindcast = hindcast, obs = obs, prThreshold = prThreshold)
    hindcast <- preprocessHindcast_res[[1]]
    minHindcastPreci <- preprocessHindcast_res[[2]]
  }
  
  # default is the simplest method in biascorrection, just do simple addition and subtraction.
  if (method == 'delta') {
    if (length(frc) != length(obs)) stop('This method needs frc data have the same length as obs data.')
    # comes from downscaleR biascorrection method
    frcMean <- mean(frc, na.rm = TRUE)
    hindcastMean <- mean(hindcast, na.rm = TRUE)
    frc <- obs - hindcastMean + frcMean
    
  } else if (method == 'scaling') {
    obsMean <- mean(obs, na.rm = TRUE)
    hindcastMean <- mean(hindcast, na.rm = TRUE)
    
    if (scaleType == 'multi') {
      frc <- frc / hindcastMean * obsMean
      
    } else if (scaleType == 'add') {
      frc <- frc - hindcastMean + obsMean
    }
    
    
  } else if (method == 'eqm') {
    if (preci == FALSE) {
      frc <- biasCorrect_core_eqm_nonPreci(frc, hindcast, obs, extrapolate, prThreshold)
    } else {
      frc <- biasCorrect_core_eqm_preci(frc, hindcast, obs, minHindcastPreci, extrapolate,
                                        prThreshold)
    }
    
  } else if (method == 'gqm') {
    if (preci == FALSE) stop ('gqm method only applys to precipitation, please set preci = T')
    frc <- biasCorrect_core_gqm(frc, hindcast, obs, prThreshold, minHindcastPreci)
  }
  
  
  return(frc)
}


#' @importFrom MASS fitdistr
#' @importFrom stats rgamma
preprocessHindcast <- function(hindcast, obs, prThreshold) {
  lowerIndex <- length(which(obs < prThreshold))
  
  # In the original function, this minHindcastPreci is Pth[,i,j] in downscaleR, and it is originally
  # set to NA, which is not so appropriate for all the precipitations.
  # In the original function, there are only two conditions, 1. all the obs less than threshold
  # 2. there are some obs less than threshold. 
  # While, if we set threshold to 0, there could be a 3rd condition, all the obs no less than threshold.
  # Here I set this situation, firstly set minHindcastPreci to the min of the hindcast. Because in future
  # use, 'eqm' method is going to use this value.
  
  # The problem above has been solved.
  
  
  if (lowerIndex >= 0 & lowerIndex < length(obs)) {
    index <- sort(hindcast, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
    hindcast_sorted <- sort(hindcast, decreasing = FALSE, na.last = NA)
    # minHindcastPreci is the min preci over threshold FOR ***HINDCAST***
    # But use obs to get the lowerIndex, so obs_sorted[lowerIndex + 1] > prThreshold, but
    # hindcast_sorted[lowerIndex + 1] may greater than or smaller than ptThreshold
    
    
    # It would be better to understand if you draw two lines: hindcast_sorted and obs_sorted
    # with y = prThreshold, you will find the difference of the two.
    
    # In principle, the value under the threshold needs to be replaced by some other reasonable value.
    # simplest way 
    minHindcastPreci <- hindcast_sorted[lowerIndex + 1]
    
    # Also here if minHindcastPreci is 0 and prThreshold is 0, will cause problem, bettter set 
    # I set it prThreshold != 0 
    if (minHindcastPreci <= prThreshold & prThreshold != 0) {
      obs_sorted <- sort(obs, decreasing = FALSE, na.last = NA)
      
      # higherIndex is based on hindcast
      higherIndex <- which(hindcast_sorted > prThreshold & !is.na(hindcast_sorted))
      
      if (length(higherIndex) == 0) {
        higherIndex <- max(which(!is.na(hindcast_sorted)))
        higherIndex <- min(length(obs_sorted), higherIndex)
      } else {
        higherIndex <- min(higherIndex)
      }
      # here I don't know why choose 6.
      # Written # [Shape parameter Scale parameter] in original package
      # according to the reference and gamma distribution, at least 6 values needed to fit gamma
      # distribution.
      if (length(unique(obs_sorted[(lowerIndex + 1):higherIndex])) < 6) {
        hindcast_sorted[(lowerIndex + 1):higherIndex] <- mean(obs_sorted[(lowerIndex + 1):higherIndex], 
                                                              na.rm = TRUE)
      } else {
        obsGamma <- fitdistr(obs_sorted[(lowerIndex + 1):higherIndex], "gamma")
        
        # this is to replace the original hindcast value between lowerIndex and higherIndex with 
        # some value taken from gamma distribution just generated.
        hindcast_sorted[(lowerIndex + 1):higherIndex] <- rgamma(higherIndex - lowerIndex, obsGamma$estimate[1], 
                                                                rate = obsGamma$estimate[2])
      }
      hindcast_sorted <- sort(hindcast_sorted, decreasing = FALSE, na.last = NA)
      
    } 
    minIndex <- min(lowerIndex, length(hindcast))
    hindcast_sorted[1:minIndex] <- 0
    hindcast[index] <- hindcast_sorted
    
  } else if (lowerIndex == length(obs)) {
    
    index <- sort(hindcast, decreasing = FALSE, na.last = NA, index.return = TRUE)$ix
    hindcast_sorted <- sort(hindcast, decreasing = FALSE, na.last = NA)
    minHindcastPreci <- hindcast_sorted[lowerIndex]
    
    # here is to compare with hindcast, not obs
    minIndex <- min(lowerIndex, length(hindcast))
    hindcast_sorted[1:minIndex] <- 0
    hindcast[index] <- hindcast_sorted
    
  }
  return(list(hindcast, minHindcastPreci))
}

biasCorrect_core_eqm_nonPreci <- function(frc, hindcast, obs, extrapolate, prThreshold) {
  ecdfHindcast <- ecdf(hindcast)
  
  if (extrapolate == 'constant') {
    higherIndex <- which(frc > max(hindcast, na.rm = TRUE))
    lowerIndex <- which(frc < min(hindcast, na.rm = TRUE))
    
    extrapolateIndex <- c(higherIndex, lowerIndex)
    non_extrapolateIndex <- setdiff(1:length(frc), extrapolateIndex)
    
    # In case extrapolateIndex is of length zero, than extrapolate cannot be used afterwards
    # So use setdiff(1:length(sim), extrapolateIndex), if extrapolateIndex == 0, than it will
    # return 1:length(sim)
    
    if (length(higherIndex) > 0) {
      maxHindcast <- max(hindcast, na.rm = TRUE)
      dif <- maxHindcast - max(obs, na.rm = TRUE)
      frc[higherIndex] <- frc[higherIndex] - dif
    }
    
    if (length(lowerIndex) > 0) {
      minHindcast <- min(hindcast, na.rm = TRUE)
      dif <- minHindcast - min(obs, nna.rm = TRUE)
      frc[lowerIndex] <- frc[lowerIndex] - dif
    }
    
    frc[non_extrapolateIndex] <- quantile(obs, probs = ecdfHindcast(frc[non_extrapolateIndex]), 
                                          na.rm = TRUE, type = 4)
  } else {
    frc <- quantile(obs, probs = ecdfHindcast(frc), na.rm = TRUE, type = 4)
  }
  return(frc)
}

biasCorrect_core_eqm_preci <- function(frc, hindcast, obs, minHindcastPreci, extrapolate, 
                                       prThreshold) {
  
  # Most of time this condition seems useless because minHindcastPreci comes from hindcast, so there will be
  # always hindcast > minHindcastPreci exists.
  # Unless one condition that minHindcastPreci is the max in the hindcast, than on hindcast > minHindcastPreci
  if (length(which(hindcast > minHindcastPreci)) > 0) {
    
    ecdfHindcast <- ecdf(hindcast[hindcast > minHindcastPreci])
    
    noRain <- which(frc <= minHindcastPreci & !is.na(frc))
    rain <- which(frc > minHindcastPreci & !is.na(frc))
    
    # drizzle is to see whether there are some precipitation between the min frc (over threshold) and 
    # min hindcast (over threshold).
    drizzle <- which(frc > minHindcastPreci & frc <= min(hindcast[hindcast > minHindcastPreci], na.rm = TRUE) 
                     & !is.na(frc))
    
    if (length(rain) > 0) {
      ecdfFrc <- ecdf(frc[rain])
      
      if (extrapolate == 'constant') {
        
        # This higher and lower index mean the extrapolation part
        higherIndex <- which(frc[rain] > max(hindcast, na.rm = TRUE))
        lowerIndex <- which(frc[rain] < min(hindcast, na.rm = TRUE))
        
        extrapolateIndex <- c(higherIndex, lowerIndex)
        non_extrapolateIndex <- setdiff(1:length(rain), extrapolateIndex)
        
        if (length(higherIndex) > 0) {
          maxHindcast <- max(hindcast, na.rm = TRUE)
          dif <- maxHindcast - max(obs, na.rm = TRUE)
          frc[rain[higherIndex]] <- frc[higherIndex] - dif
        }
        
        if (length(lowerIndex) > 0) {
          minHindcast <- min(hindcast, na.rm = TRUE)
          dif <- minHindcast - min(obs, nna.rm = TRUE)
          frc[rain[lowerIndex]] <- frc[lowerIndex] - dif
        }
        
        # Here the original function doesn't accout for the situation that extraploateIndex is 0
        # if it is 0, rain[-extraploateIndex] would be nothing
        
        # Above has been solved by using setdiff.
        frc[rain[non_extrapolateIndex]] <- quantile(obs[which(obs > prThreshold & !is.na(obs))], 
                                                    probs = ecdfHindcast(frc[rain[non_extrapolateIndex]]), 
                                                    na.rm = TRUE, type = 4)
      } else {
        
        frc[rain] <- quantile(obs[which(obs > prThreshold & !is.na(obs))], 
                              probs = ecdfHindcast(frc[rain]), na.rm = TRUE, type = 4)
      }
    }
    if (length(drizzle) > 0){
      
      # drizzle part is a seperate part. it use the ecdf of frc (larger than minHindcastPreci) to 
      # biascorrect the original drizzle part
      frc[drizzle] <- quantile(frc[which(frc > min(hindcast[which(hindcast > minHindcastPreci)], na.rm = TRUE) & 
                                           !is.na(frc))], probs = ecdfFrc(frc[drizzle]), na.rm = TRUE, 
                               type = 4)
    }
    
    frc[noRain] <- 0
    
  } else {
    # in this condition minHindcastPreci is the max of hindcast, so all hindcast <= minHindcastPreci
    # And frc distribution is used then.
    noRain <- which(frc <= minHindcastPreci & !is.na(frc))
    rain <- which(frc > minHindcastPreci & !is.na(frc))
    
    if (length(rain) > 0) {
      ecdfFrc <- ecdf(frc[rain])
      frc[rain] <- quantile(obs[which(obs > prThreshold & !is.na(obs))], probs = ecdfFrc(frc[rain]), 
                            na.rm = TRUE, type = 4)
    }
    frc[noRain]<-0
  }
  return(frc)
}

biasCorrect_core_gqm <- function(frc, hindcast, obs, prThreshold, minHindcastPreci) {
  if (any(obs > prThreshold)) {
    
    ind <- which(obs > prThreshold & !is.na(obs))
    obsGamma <- fitdistr(obs[ind],"gamma")
    ind <- which(hindcast > 0 & !is.na(hindcast))
    hindcastGamma <- fitdistr(hindcast[ind],"gamma")
    rain <- which(frc > minHindcastPreci & !is.na(frc))
    noRain <- which(frc <= minHindcastPreci & !is.na(frc))
    
    probF <- pgamma(frc[rain], hindcastGamma$estimate[1], rate = hindcastGamma$estimate[2])
    frc[rain] <- qgamma(probF,obsGamma$estimate[1], rate = obsGamma$estimate[2])
    frc[noRain] <- 0
  } else {
    warning('All the observations of this cell(station) are lower than the threshold, 
            no bias correction applied.')
  }
  return(frc)
  }
