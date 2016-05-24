



#' Get bias factor for multi/operational/real time bias correction.
#' 
#' When you do multi/operational/real time bias correction. It's too expensive
#' to input hindcast and obs every time. Especially when you have a long period of hindcast
#' and obs, but only a short period of frc, it's too unecessary to read and compute hindcast
#' and obs everytime. Therefore, biasFactor is designed. Using \code{getBiasFactor}, you can
#' get the biasFactor with hindcast and observation, then you can use \code{applyBiasFactor} to 
#' apply the biasFactor to different forecasts. 
#' 
#' @param hindcast a hyfo grid data output or a dataframe(time series) consists of Date column and one or more value columns, 
#' representing the hindcast data. This data will be used in the calibration of the forecast, so it's better to have the same date period as
#' observation data. Check details for more information.
#' @param obs a hyfo grid data output or a dataframe (time series) consists of Date column and one or more value columns, 
#' representing the observation data.
#' @param method bias correct method, including 'delta', 'scaling'...,default method is 'scaling'.
#' @param scaleType only when the method "scaling" is chosen, scaleType will be available. Two different types
#' of scaling method, 'add' and 'multi', which means additive and multiplicative scaling method, default is 'multi'. More info check 
#' details.
#' @param preci If the precipitation is biascorrected, then you have to assign \code{preci = TRUE}. Since for
#' precipitation, some biascorrect methods may not apply to, or some methods are specially for precipitation. 
#' Default is FALSE, refer to details.
#' @param prThreshold The minimum value that is considered as a non-zero precipitation. Default to 1 (assuming mm).
#' @param extrapolate When use 'eqm' method, and 'no' is set, modified frc is bounded by the range of obs.
#' If 'constant' is set, modified frc is not bounded by the range of obs. Default is 'no'.
#' 
#' @seealso \code{\link{biasCorrect}} for method used in bias correction.
#' \code{\link{applyBiasFactor}}, for the second part.
#' 
#' @details 
#' 
#' Information about the method and how biasCorrect works can be found in \code{\link{biasCorrect}}
#' 
#' \strong{why use biasFactor}
#' 
#' As for forecasting, for daily data, there is usually no need to have
#' different bias factor every different day. You can calculate one bisa factor using a long
#' period of hindcast and obs, and apply that factor to different frc.
#' 
#' For example,
#' 
#' You have 10 years of hindcast and observation. you want to do bias correction for some 
#' forecasting product, e.g. system 4. For system 4, each month, you will get a new forecast
#' about the future 6 months. So if you want to do the real time bias correction, you have to
#' take the 10 years of hindcast and observation data with you, and run \code{biasCorrect} every
#' time you get a new forecast. That's too expensive.
#' 
#' For some practical use in forecasting, there isn't a so high demand for accuracy. E.g.,
#' Maybe for February and March, you can use the same biasFactor, no need to do the computation 
#' again. 
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
#' # Then we will use nc data as forecasting data, and use itself as hindcast data,
#' # use tgridData as observation.
#' 
#' biasFactor <- getBiasFactor(nc, tgridData)
#' newFrc <- applyBiasFactor(nc, biasFactor)
#'    
#' biasFactor <- getBiasFactor(nc, tgridData, method = 'eqm', extrapolate = 'constant',
#' preci = TRUE)
#' # This method needs obs input.
#' newFrc <- applyBiasFactor(nc, biasFactor, obs = tgridData)
#' 
#' biasFactor <- getBiasFactor(nc, tgridData, method = 'gqm', preci = TRUE)
#' newFrc <- applyBiasFactor(nc, biasFactor) 
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
#' # default method is scaling
#' biasFactor <- getBiasFactor(hindcast, obs)
#' frc_new <- applyBiasFactor(frc, biasFactor)
#' 
#' # for precipitation data, extra process needs to be executed, so you have to tell
#' # the program to it is a precipitation data.
#' 
#' biasFactor <- getBiasFactor(hindcast, obs, preci = TRUE)
#' frc_new1 <- applyBiasFactor(frc, biasFactor)
#' 
#' # You can use other methods to biascorrect, e.g. delta method. 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'delta')
#' # delta method needs obs input.
#' frc_new2 <- applyBiasFactor(frc, biasFactor, obs = obs)
#' 
#' # 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'eqm', preci = TRUE)
#' # eqm needs obs input
#' frc_new3 <- applyBiasFactor(frc, biasFactor, obs = obs)
#' 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'gqm', preci = TRUE)
#' frc_new4 <- applyBiasFactor(frc, biasFactor)
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
#' 
#' @importFrom methods setMethod
#' @export
#' 
#' 
# debug by trace("getBiasFactor", browser, exit=browser, signature = c("list", "list"))
setGeneric('getBiasFactor', function(hindcast, obs, method = 'scaling', scaleType = 'multi', 
                                     preci = FALSE, prThreshold = 0, extrapolate = 'no') {
  standardGeneric('getBiasFactor')
})

#' @describeIn getBiasFactor
setMethod('getBiasFactor', signature('data.frame', 'data.frame'), 
          function(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
            result <- getBiasFactor.TS(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate)
            return(result)
          })


# This is for the grid file from downscaleR
#' @describeIn getBiasFactor
#' @importFrom methods new
setMethod('getBiasFactor', signature('list', 'list'), 
          function(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
            result <- getBiasFactor.list(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate)
            return(result)
            })




#' Apply bias factor to different forecasts for multi/operational/real time bias correction.
#' 
#' When you do multi/operational/real time bias correction. It's too expensive
#' to input hindcast and obs every time. Especially when you have a long period of hindcast
#' and obs, but only a short period of frc, it's too unecessary to read and compute hindcast
#' and obs everytime. Therefore, biasFactor is designed. Using \code{getBiasFactor}, you can
#' get the biasFactor with hindcast and observation, then you can use \code{applyBiasFactor} to 
#' apply the biasFactor to different forecasts. 
#' 
#' @param frc a hyfo grid data output or a dataframe(time series) consists of Date column and one or more value columns, 
#' representing the frc data. Check details for more information.
#' @param biasFactor a file containing all the information of the calibration, will be
#' applied to different forecasts.
#' @param obs for some methods, observation input is necessary. obs is a hyfo grid data output or a dataframe (time series) consists of Date column and one or more value columns, 
#' representing the observation data. Default value is NULL.
#' @seealso \code{\link{biasCorrect}} for method used in bias correction. 
#' \code{\link{getBiasFactor}}, for the first part.
#' 
#' @details 
#' 
#' Information about the method and how biasCorrect works can be found in \code{\link{biasCorrect}}
#' 
#' \strong{why use biasFactor}
#' 
#' As for forecasting, for daily data, there is usually no need to have
#' different bias factor every different day. You can calculate one bisa factor using a long
#' period of hindcast and obs, and apply that factor to different frc.
#' 
#' For example,
#' 
#' You have 10 years of hindcast and observation. you want to do bias correction for some 
#' forecasting product, e.g. system 4. For system 4, each month, you will get a new forecast
#' about the future 6 months. So if you want to do the real time bias correction, you have to
#' take the 10 years of hindcast and observation data with you, and run \code{biasCorrect} every
#' time you get a new forecast. That's too expensive.
#' 
#' For some practical use in forecasting, there isn't a so high demand for accuracy. E.g.,
#' Maybe for February and March, you can use the same biasFactor, no need to do the computation 
#' again. 
#' 
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
#' #' # Since the example data, has some NA values, the process will include some warning #message, 
#' # which can be ignored in this case.
#' 
#' 
#' 
#' # Then we will use nc data as forecasting data, and use itself as hindcast data,
#' # use tgridData as observation.
#' 
#' biasFactor <- getBiasFactor(nc, tgridData)
#' newFrc <- applyBiasFactor(nc, biasFactor)
#'    
#' biasFactor <- getBiasFactor(nc, tgridData, method = 'eqm', extrapolate = 'constant',
#' preci = TRUE)
#' # This method needs obs input.
#' newFrc <- applyBiasFactor(nc, biasFactor, obs = tgridData)
#' 
#' biasFactor <- getBiasFactor(nc, tgridData, method = 'gqm', preci = TRUE)
#' newFrc <- applyBiasFactor(nc, biasFactor) 
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
#' # default method is scaling
#' biasFactor <- getBiasFactor(hindcast, obs)
#' frc_new <- applyBiasFactor(frc, biasFactor)
#' 
#' # for precipitation data, extra process needs to be executed, so you have to tell
#' # the program to it is a precipitation data.
#' 
#' biasFactor <- getBiasFactor(hindcast, obs, preci = TRUE)
#' frc_new1 <- applyBiasFactor(frc, biasFactor)
#' 
#' # You can use other methods to biascorrect, e.g. delta method. 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'delta')
#' # delta method needs obs input.
#' frc_new2 <- applyBiasFactor(frc, biasFactor, obs = obs)
#' 
#' # 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'eqm', preci = TRUE)
#' # eqm needs obs input
#' frc_new3 <- applyBiasFactor(frc, biasFactor, obs = obs)
#' 
#' biasFactor <- getBiasFactor(hindcast, obs, method = 'gqm', preci = TRUE)
#' frc_new4 <- applyBiasFactor(frc, biasFactor)
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
#' 
#' @export
setGeneric('applyBiasFactor', function(frc, biasFactor, obs = NULL) {
  standardGeneric('applyBiasFactor')
})

#' @describeIn applyBiasFactor
#' @importFrom methods setMethod
setMethod('applyBiasFactor', signature('data.frame', 'biasFactor'), 
          function(frc, biasFactor, obs) {
            result <- applyBiasFactor.TS(frc, biasFactor, obs)
            return(result)
          })
           
#' @describeIn applyBiasFactor
#' @importFrom methods setMethod
setMethod('applyBiasFactor', signature('list', 'biasFactor.hyfo'), 
          function(frc, biasFactor, obs) {
            result <- applyBiasFactor.list(frc, biasFactor, obs)
            return(result)
          })


### generic functions
getBiasFactor.TS <- function(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
  
  if (!grepl('-|/', obs[1, 1]) | !grepl('-|/', hindcast[1, 1])) {
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
  n <- ncol(hindcast)
  
  # For every column, it's biascorrected respectively.
  biasFactor <- lapply(2:n, function(x) getBiasFactor_core(hindcast[, x], obs[, 2], method = method,
                                                           scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                           extrapolate = extrapolate))
  if (n - 1 > 1) {
    biasFactor_all <- new('biasFactor.multiMember', biasFactor = biasFactor, memberDim = n - 1,
                          method = method, preci = preci, prThreshold = prThreshold, scaleType = scaleType, 
                          extrapolate = extrapolate)
    
  } else {
    biasFactor_all <- new('biasFactor', biasFactor = biasFactor, method = method, 
                          preci = preci, prThreshold = prThreshold, scaleType = scaleType, 
                          extrapolate = extrapolate)
  }
  
  return(biasFactor_all)
}

getBiasFactor.list <- function(hindcast, obs, method, scaleType, preci, prThreshold, extrapolate) {
  
  ## Check if the data is a hyfo grid data.
  checkHyfo(hindcast, obs)
  
  hindcastData <- hindcast$Data
  obsData <- obs$Data
  
  ## save frc dimension order, at last, set the dimension back to original dimension
  hindcastDim <- attributes(hindcastData)$dimensions
  
  ## ajust the dimension into general dimension order.
  obsData <- adjustDim(obsData, ref = c('lon', 'lat', 'time'))
  
  ## CheckDimLength, check if all the input dataset has different dimension length
  # i.e. if they all have the same lon and lat number.
  checkDimLength(hindcastData, obsData, dim = c('lon', 'lat'))
  
  
  # Now real bias correction is executed.
  
  memberIndex <- grepAndMatch('member', attributes(hindcastData)$dimensions)
  
  # For dataset that has a member part 
  if (!is.na(memberIndex)) {
    
    hindcastData <- adjustDim(hindcastData, ref = c('lon', 'lat', 'time', 'member'))
    
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
    
    biasFactor_all <- vector(mode = "list", length = dim(hindcastData)[4])
    for (member in 1:dim(hindcastData)[4]) {
      biasFactor_all[[member]] <- vector(mode = 'list', length = dim(hindcastData)[1])
      for (lon in 1:dim(hindcastData)[1]) {
        biasFactor_all[[member]][[lon]] <- vector(mode = 'list', length = dim(hindcastData)[2])
        for (lat in 1:dim(hindcastData)[2]) {
          biasFactor_all[[member]][[lon]][[lat]] <- getBiasFactor_core(hindcastData[lon, lat,, member], obsData[lon, lat,], method = method,
                                                                       scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                                       extrapolate = extrapolate)
        }
      }
    }
    
    biasFactor <- new('biasFactor.hyfo', biasFactor = biasFactor_all, method = method, preci = preci,
                      prThreshold = prThreshold, scaleType = scaleType, extrapolate = extrapolate, 
                      lonLatDim = calcuDim(hindcastData, dim = c('lon', 'lat')),
                      memberDim = calcuDim(hindcastData, dim = 'member'))
  } else {
    
    hindcastData <- adjustDim(hindcastData, ref = c('lon', 'lat', 'time'))
    
    biasFactor_all <- vector(mode = 'list', length = dim(hindcastData)[1])
    for (lon in 1:dim(hindcastData)[1]) {
      biasFactor_all[[lon]] <- vector(mode = 'list', length = dim(hindcastData)[2]) 
      for (lat in 1:dim(hindcastData)[2]) {
        biasFactor_all[[lon]][[lat]] <- getBiasFactor_core(hindcastData[lon, lat,], obsData[lon, lat,], method = method,
                                                           scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                           extrapolate = extrapolate)
      }
    }
    biasFactor <- new('biasFactor.hyfo', biasFactor = biasFactor_all, method = method, preci = preci,
                      prThreshold = prThreshold, scaleType = scaleType, extrapolate = extrapolate, 
                      lonLatDim = calcuDim(hindcastData, dim = c('lon', 'lat')))
    
  }
  
  return(biasFactor)
}

applyBiasFactor.TS <- function(frc, biasFactor, obs) {
  method <- biasFactor@method
  preci <- biasFactor@preci
  prThreshold <- biasFactor@prThreshold
  scaleType <- biasFactor@scaleType
  extrapolate <- biasFactor@extrapolate
  memberDim <- biasFactor@memberDim
  biasFactor <- biasFactor@biasFactor
  
  
  # First check if the first column is Date
  if (!grepl('-|/', frc[1, 1])) {
    stop('First column is not date or Wrong Date formate, check the format in ?as.Date{base} 
         and use as.Date to convert.If your input is a hyfo dataset, put input = "hyfo" as an
         argument, check help for more info.')
  }
  # change to date type is easier, but in case in future the flood part is added, Date type doesn't have
  # hour, min and sec, so, it's better to convert it into POSIxlt.
  
  # In this case more than one value columns exist in the dataset, both frc and hindcast.
  
  n <- ncol(frc)
  if (n-1 != memberDim) stop('frc and biasFactor have different members.')
  
  
  # For every column, it's biascorrected respectively.
  frc_data <- lapply(2:n, function(x) applyBiasFactor_core(frc[, x], biasFactor = biasFactor[[x - 1]], method = method,
                                                           scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                           extrapolate = extrapolate, obs = obs[, 2]))
  frc_data <- do.call('cbind', frc_data)
  rownames(frc_data) <- NULL
  
  names <- colnames(frc)
  frc_new <- data.frame(frc[, 1], frc_data)
  colnames(frc_new) <- names
  
  return(frc_new)
  
}

applyBiasFactor.list <- function(frc, biasFactor, obs) {
  method <- biasFactor@method
  preci <- biasFactor@preci
  prThreshold <- biasFactor@prThreshold
  scaleType <- biasFactor@scaleType
  extrapolate <- biasFactor@extrapolate
  lonLatDim <- biasFactor@lonLatDim
  memberDim <- biasFactor@memberDim
  biasFactor <- biasFactor@biasFactor
  
  ## Check if the data is a hyfo grid data.
  checkHyfo(frc)
  
  
  obsData <- obs$Data
  frcData <- frc$Data
  
  ## save frc dimension order, at last, set the dimension back to original dimension
  frcDim <- attributes(frcData)$dimensions
  
  ## ajust the dimension into general dimension order.
  obsData <- adjustDim(obsData, ref = c('lon', 'lat', 'time'))
  
  ## CheckDimLength, check if all the input dataset has different dimension length
  # i.e. if they all have the same lon and lat number.
  if (!identical(calcuDim(frcData, dim = c('lon', 'lat')), lonLatDim)) {
    stop('frc data has different lon and lat from hindcast data.')
  }
  
  
  # Now real bias correction is executed.
  
  memberIndex <- grepAndMatch('member', attributes(frcData)$dimensions)
  
  # For dataset that has a member part 
  if (!is.na(memberIndex)) {
    # check if frcData and hindcastData has the same dimension and length.
    if (calcuDim(frcData, dim = 'member') != memberDim) {
      stop('frc data has different member number from hindcast.')
    } 
    
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
          frcData[lon, lat,, member] <- applyBiasFactor_core(frcData[lon, lat,,member], biasFactor = biasFactor[[member]][[lon]][[lat]], method = method,
                                                             scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                             extrapolate = extrapolate, obs = obsData[lon, lat,])
        }
      }
    }
  } else {
    frcData <- adjustDim(frcData, ref = c('lon', 'lat', 'time'))
    for (lon in 1:dim(frcData)[1]) {
      for (lat in 1:dim(frcData)[2]) {
        frcData[lon, lat,] <- applyBiasFactor_core(frcData[lon, lat,], biasFactor = biasFactor[[lon]][[lat]], method = method,
                                                   scaleType = scaleType, preci = preci, prThreshold = prThreshold, 
                                                   extrapolate = extrapolate, obs = obsData[lon, lat,])
      }
    }
  }
  
  frcData <- adjustDim(frcData, ref = frcDim)
  frc$Data <- frcData
  frc$biasCorrected_by <- method
  frc_new <- frc
  
  return(frc_new)
}


#################
################# core functions for multi bias correction.

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
getBiasFactor_core <- function(hindcast, obs, method , scaleType, preci, prThreshold, extrapolate){
  # If the variable is precipitation, some further process needs to be added.
  # The process is taken from downscaleR, to provide a more reasonable hindcast, used in the calibration.
  
  
  # check if frc, hindcast or obs are all na values
  if (!any(!is.na(obs)) | !any(!is.na(hindcast))) {
    warning('In this cell, hindcast or obs data is missing. No biasCorrection for this cell.')
    return(NA)
  } 
  
  if (preci == TRUE) {
    preprocessHindcast_res <- preprocessHindcast(hindcast = hindcast, obs = obs, prThreshold = prThreshold)
    hindcast <- preprocessHindcast_res[[1]]
    minHindcastPreci <- preprocessHindcast_res[[2]]
  }
  
  # default is the simplest method in biascorrection, just do simple addition and subtraction.
  if (method == 'delta') {
    biasFactor <- getBiasFactor_core_delta(hindcast)
  } else if (method == 'scaling') {
    biasFactor <- getBiasFactor_core_scaling(hindcast, obs, scaleType)
  } else if (method == 'eqm') {
    # In this method, the value is bounded by the observation
    # Preci or not both have the same biasFactor
    if (preci == FALSE) {
      biasFactor <- getBiasFactor_core_eqm_nonPreci(hindcast, obs, extrapolate)
    } else {
      biasFactor <- getBiasFactor_core_eqm_preci(hindcast, obs, minHindcastPreci, extrapolate, prThreshold)
    }
    
    
  } else if (method == 'gqm') {
    if (preci == FALSE) stop ('gqm method only applys to precipitation, please set preci = T')
    biasFactor <- getBiasFactor_core_gqm(hindcast, obs, prThreshold, minHindcastPreci)
  }
  
  if (preci == TRUE) biasFactor$minHindcastPreci <- minHindcastPreci
  
  return(biasFactor)
}


applyBiasFactor_core <- function(frc, biasFactor, method, preci, prThreshold, scaleType,
                                 extrapolate, obs = NULL) {
  
  if (!any(!is.na(biasFactor))) {
    warning('In this cell, biasFactor is missing.No biasCorrection for this cell.')
    # here return NA or return the unprocessed frc, both are OK. But return NA is more
    # obvious for user.
    return(NA)
  }
  
  if (method == 'delta') {
    if (is.null(obs)) stop('This method needs obs input.')
    if (length(frc) != length(obs)) stop('This method needs frc data have the same length as obs data.')
    frc <- applyBiasFactor_core_delta(frc = frc, biasFactor = biasFactor, obs = obs)
  } else if (method == 'scaling') {
    frc <- applyBiasFactor_core_scaling(frc = frc, biasFactor = biasFactor, scaleType = scaleType)
  } else if (method == 'eqm') {
    if (is.null(obs)) stop('This method needs obs input.')
    if (preci == FALSE) {
      frc <- applyBiasFactor_core_eqm_nonPreci(frc = frc, biasFactor = biasFactor, extrapolate = extrapolate, 
                                               obs = obs)
    } else {
      frc <- applyBiasFactor_core_eqm_preci(frc = frc, biasFactor = biasFactor, extrapolate = extrapolate, 
                                            prThreshold = prThreshold, obs = obs)
    }
  } else if (method == 'gqm') {
    frc <- applyBiasFactor_core_gqm(frc = frc, biasFactor = biasFactor)
  }
  
  return(frc)
}


getBiasFactor_core_delta <- function(hindcast) {
  biasFactor <- list()
  biasFactor$hindcastMean <- mean(hindcast, na.rm = TRUE)
  return(biasFactor)
}
applyBiasFactor_core_delta <- function(frc, biasFactor, obs) {
  hindcastMean <- biasFactor$hindcastMean
  frcMean <- mean(frc, na.rm = TRUE)
  return(obs - hindcastMean + frcMean)
}

getBiasFactor_core_scaling <- function(hindcast, obs, scaleType) {
  biasFactor <- list()
  
  hindcastMean <- mean(hindcast, na.rm = TRUE)
  obsMean <- mean(obs, na.rm = TRUE)
  
  if (scaleType == 'multi') {
    biasFactor$scale <- obsMean / hindcastMean
    
  } else if (scaleType == 'add') {
    biasFactor$scale <- obsMean - hindcastMean
  }
  
  return(biasFactor)
}

applyBiasFactor_core_scaling <- function(frc, biasFactor, scaleType) {
  
  if (scaleType == 'multi') {
    frc <- frc * biasFactor$scale
    
  } else if (scaleType == 'add') {
    frc <- frc + biasFactor$scale
  }
  return(frc)
}

getBiasFactor_core_eqm_nonPreci <- function(hindcast, obs, extrapolate) {
  
  biasFactor <- list()
  biasFactor$ecdfHindcast <- ecdf(hindcast)
  
  if (extrapolate == 'constant') {
    biasFactor$maxHindcast <- max(hindcast, na.rm = TRUE)
    biasFactor$minHindcast <- min(hindcast, na.rm = TRUE)
    biasFactor$higherIndex_dif <- biasFactor$maxHindcast - max(obs, na.rm = TRUE)
    biasFactor$lowerIndex_dif <- biasFactor$minHindcast - min(obs, na.rm = TRUE)
  }
  
  return(biasFactor)
}

getBiasFactor_core_eqm_preci <- function(hindcast, obs, minHindcastPreci, extrapolate,
                                         prThreshold) {
  
  biasFactor <- list()
  biasFactor$ecdfHindcast <- ecdf(hindcast[hindcast > minHindcastPreci])
  
  if (extrapolate == 'constant') {
    biasFactor$maxHindcast <- max(hindcast, na.rm = TRUE)
    biasFactor$minHindcast <- min(hindcast, na.rm = TRUE)
    biasFactor$higherIndex_dif <- biasFactor$maxHindcast - max(obs, na.rm = TRUE)
    biasFactor$lowerIndex_dif <- biasFactor$minHindcast - min(obs, nna.rm = TRUE)
  }
  biasFactor$availableHindcastLength <- length(which(hindcast > minHindcastPreci))
  
  # drizzle parameter 1
  biasFactor$drizzleP1 <- min(hindcast[hindcast > minHindcastPreci], na.rm = TRUE)
  # biasFactor$prThreshold <- prThreshold
  return(biasFactor)
}

applyBiasFactor_core_eqm_nonPreci <- function(frc, biasFactor, extrapolate, obs) {
  ecdfHindcast <- biasFactor$ecdfHindcast
  
  if (extrapolate == 'constant') {
    higherIndex <- which(frc > biasFactor$maxHindcast)
    lowerIndex <- which(frc < biasFactor$minHindcast)
    
    extrapolateIndex <- c(higherIndex, lowerIndex)
    non_extrapolateIndex <- setdiff(1:length(frc), extrapolateIndex)
    
    # In case extrapolateIndex is of length zero, than extrapolate cannot be used afterwards
    # So use setdiff(1:length(sim), extrapolateIndex), if extrapolateIndex == 0, than it will
    # return 1:length(sim)
    
    if (length(higherIndex) > 0) {
      
      frc[higherIndex] <- frc[higherIndex] - biasFactor$higherIndex_dif
    }
    
    if (length(lowerIndex) > 0) {
      
      frc[lowerIndex] <- frc[lowerIndex] - biasFactor$lowerIndex_dif
    }
    
    frc[non_extrapolateIndex] <- quantile(obs, probs = ecdfHindcast(frc[non_extrapolateIndex]), 
                                          na.rm = TRUE, type = 4)
  } else {
    frc <- quantile(obs, probs = ecdfHindcast(frc), na.rm = TRUE, type = 4)
  }
  return(frc)
}

#' @importFrom stats quantile
applyBiasFactor_core_eqm_preci <- function(frc, biasFactor, extrapolate, prThreshold, obs) {
  
  # Most of time this condition seems useless because minHindcastPreci comes from hindcast, so there will be
  # always hindcast > minHindcastPreci exists.
  # Unless one condition that minHindcastPreci is the max in the hindcast, than on hindcast > minHindcastPreci
  if (biasFactor$availableHindcastLength > 0) {
    
    ecdfHindcast <- biasFactor$ecdfHindcast
    
    noRain <- which(frc <= biasFactor$minHindcastPreci & !is.na(frc))
    rain <- which(frc > biasFactor$minHindcastPreci & !is.na(frc))
    
    # drizzle is to see whether there are some precipitation between the min frc (over threshold) and 
    # min hindcast (over threshold).
    drizzle <- which(frc > biasFactor$minHindcastPreci & frc <= biasFactor$drizzleP1 & !is.na(frc))
    
    if (length(rain) > 0) {
      ecdfFrc <- ecdf(frc[rain])
      
      if (extrapolate == 'constant') {
        
        # This higher and lower index mean the extrapolation part
        higherIndex <- which(frc[rain] > biasFactor$maxHindcast)
        lowerIndex <- which(frc[rain] < biasFactor$minHindcast)
        
        extrapolateIndex <- c(higherIndex, lowerIndex)
        non_extrapolateIndex <- setdiff(1:length(rain), extrapolateIndex)
        
        if (length(higherIndex) > 0) {
          frc[rain[higherIndex]] <- frc[higherIndex] - biasFactor$higherIndex_dif
        }
        
        if (length(lowerIndex) > 0) {
          frc[rain[lowerIndex]] <- frc[lowerIndex] - biasFactor$lowerIndex_dif
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
      frc[drizzle] <- quantile(frc[which(frc > biasFactor$drizzleP1 & !is.na(frc))], 
                               probs = ecdfFrc(frc[drizzle]), na.rm = TRUE, 
                               type = 4)
    }
    
    frc[noRain] <- 0
    
  } else {
    # in this condition minHindcastPreci is the max of hindcast, so all hindcast <= minHindcastPreci
    # And frc distribution is used then.
    noRain <- which(frc <= biasFactor$minHindcastPreci & !is.na(frc))
    rain <- which(frc > biasFactor$minHindcastPreci & !is.na(frc))
    
    if (length(rain) > 0) {
      ecdfFrc <- ecdf(frc[rain])
      frc[rain] <- quantile(obs[which(obs > prThreshold & !is.na(obs))], probs = ecdfFrc(frc[rain]), 
                            na.rm = TRUE, type = 4)
    }
    frc[noRain]<-0
  }
  return(frc)
}

#' @importFrom MASS fitdistr
getBiasFactor_core_gqm <- function(hindcast, obs, prThreshold, minHindcastPreci) {
  if (any(obs > prThreshold)) {
    biasFactor <- list()
    ind <- which(obs > prThreshold & !is.na(obs))
    obsGamma <- fitdistr(obs[ind],"gamma")
    biasFactor$obsShape <- obsGamma$estimate[1]
    biasFactor$obsRate <- obsGamma$estimate[2]
    
    ind <- which(hindcast > 0 & !is.na(hindcast))
    hindcastGamma <- fitdistr(hindcast[ind],"gamma")
    biasFactor$hindcastShape <- hindcastGamma$estimate[1]
    biasFactor$hindcastRate <- hindcastGamma$estimate[2]
    biasFactor$minHindcastPreci <- minHindcastPreci
    
  } else {
    warning('All the observations of this cell(station) are lower than the threshold, 
            no biasFactor returned.')
    biasFactor <- NA
  }
  return(biasFactor)
}

#' @importFrom stats pgamma qgamma
applyBiasFactor_core_gqm <- function(frc, biasFactor) {
  
  rain <- which(frc > biasFactor$minHindcastPreci & !is.na(frc))
  noRain <- which(frc <= biasFactor$minHindcastPreci & !is.na(frc))
  
  probF <- pgamma(frc[rain], biasFactor$hindcastShape, rate = biasFactor$hindcastRate)
  frc[rain] <- qgamma(probF, biasFactor$obsShape, rate = biasFactor$obsRate)
  frc[noRain] <- 0
  
  return(frc)
}