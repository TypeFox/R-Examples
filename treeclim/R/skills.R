##' Evaluate reconstruction skills using split-calibration
##' 
##' This function allows to evaluate the reconstruction skills for a
##' given proxy time series in split-calibration approach.
##' @details The result of a call to \code{\link{dcc}} or
##' \code{\link{seascorr}} can be used as \code{object} for the
##' function. The required data is then taken from this object and no
##' further processing of the tree and climate variables has to be
##' done for by the user. This reflects the flow of analysis, where
##' first general climate/growth relations are explored, and then the
##' strongest ones are deployed for reconstruction purposes.
##'
##' \code{target} is an aggregation modifier (one of
##' \code{\link{.mean}}, \code{\link{.sum}}, and
##' \code{\link{.range}}). The user should be aware of the fact that
##' in case the aggregation modifier evaluates to more than one
##' variable (e.g., summer means for both temperature and
##' precipiation), a warning message is issued, and only the first
##' variable is taken into consideration for evaluating the
##' reconstruction skills. If not specified, the selection from the
##' original call to \code{\link{dcc}} is used.
##'
##' The type of regression model (ordinary least squares or
##' errors-in-variables via reduced major axis regression) can be
##' selected.
##'
##' The part of the data to be used as a calibration subset can be
##' specified in three different ways: 1) as a range of years, these
##' are then taken as calibration period; 2) as a single integer, if
##' positive, this number of observations at the recent end of the
##' data set is taken as calibration set, if negative, this number of
##' oldest observations is taken; and 3) as a character string giving
##' a percentage of values, e.g., "-40\%" would select the 40\% oldest
##' observations, while "55\%" would select the 55\% most recent ones.
##'
##' The relationship between climate and tree-ring data is evaluated
##' for the calibration period and the complete data set. Frequently
##' used verification statistics are computed: reduction of error
##' (RE), coefficient of efficiency (CE), and the Durban-Watson
##' statistic (DW) (Cook et al. 1994, Durbin and Watson, 1951).
##' @param object an object of class "tc_dcc" or "tc_seascorr"
##' @param target a treeclim selection modifier specifying the climate
##'     target to be reconstructed, see below for details
##' @param model one of "ols" or "rma"
##' @param calibration which part of the data shall be used as
##'     calibration subset? Given as either a range of years, an
##'     integer corresponding to the first or last number of
##'     observations, or a percentage as character string
##'     corresponding to the part of the data set to be used as
##'     calibration subset.
##' @param timespan timespan to be used to truncate the data
##' @return 'skills' returns an 'object' of class '"tc_skills"'.
##' 
##' An object of class '"tc_skills"' is a list containing at least the
##' following components:
##' 
##' \item{call}{the call made to function 'skills'}
##' \item{target}{the target used for reconstruction}     
##' \item{r.cal}{the coefficient of correlation for the calibration
##' timespan}    
##' \item{r.full}{the coefficient of correlation for the complete data
##' set}   
##' \item{coef.cal}{regression coefficients for the calibration model} 
##' \item{coef.full}{regression coefficients for the full model}
##' \item{p.cal}{significance for the calibration model}    
##' \item{p.full}{significance for the full model}   
##' \item{RE}{reduction of error statistic}       
##' \item{CE}{coefficient of efficiency statistic}       
##' \item{DW}{Durbin-Watson statistic}       
##' \item{cal.model}{the complete calibration model (an object of
##' class 'lmodel2')}
##' \item{full.mode}{the complete full model (an object of class 'lmodel2')}
##' @references
##' Cook E, Briffa K, Jones P (1994) Spatial regression
##' methods in dendroclimatology: A review and comparison of two
##' techniques. International Journal of Climatology, 14, 379-402.
##'
##' Durbin, J, Watson, GS (1951) Testing for serial
##' correlation in least squares regression. Biometrika 38:159-78.
##' @examples
##' \dontrun{
##' dc <- dcc(muc_fake, muc_clim, .mean(6:9, "temp") + .sum(6:9,
##' "prec"))
##' sk <- skills(dc)
##' sk
##' plot(sk)
##' }
##' @import lmtest
##' @import lmodel2
##' @export
skills <- function(object, target = NULL, model = "ols",
                  calibration = "50%", timespan = NULL) {

  if (!any(c("tc_dcc", "tc_seascorr") == class(object)[1])) {
    stop("`object` must be the result of functions `dcc` or `seascorr`.")
  }

  Method <- NULL                        # to keep R CMD check happy
  
  mf <- match.call()
  if (is.null(mf$target)) {
      mf$target <- object$call$selection
      if (is.null(mf$target)) {
        mf$target <- -6:9
      }
  }
  x_sel <- eval(mf$target)
  if (is.numeric(x_sel)) {
    x_sel <- eval(
      substitute(
        .range(.months = sel, .variables = NULL),
        list(sel = x_sel)
      )
    )
  }
    
  monthcheck <- check_months(eval(x_sel))
  if (monthcheck$check == FALSE) {
    stop("Please specify months with numbers from -12 (previous december) to 12 (current december).")
  }
  minmonth <- monthcheck$minmonth
  
  truncated_data <-
      truncate_input(object$original$tree,
                     object$original$climate, timespan,
                     minmonth, FALSE, TRUE)
  pad <- truncated_data$pad
  climate <- truncated_data$climate
  chrono <- truncated_data$chrono
    
  pmat <- make_pmat(climate, pad)
      
  X <- tc_design(eval(x_sel), pmat, check_2 = FALSE)

  all_years <- as.numeric(rownames(X$aggregate))
  m <- length(all_years)
  
  if (dim(X$aggregate)[2] > 1) {
    warning(paste0("Reconstruction skills cannot be evaluated when using more than one independent variable. We use only the first variable (by alphabet: ", X$names[1], ")."))
  }

  x <- X$aggregate[,1]

  ## split into calibration and verification periods first, deparse
  ## what is given as calibration period

  if (is.character(calibration)) {
    if (length(grep("^-{0,1}[0-9]+%$", calibration)) == 1) {
      ## calibration given as percentage
      num <- gregexpr("[0-9]+", calibration)[[1]]
      perc <- substr(calibration, num,
                     num - 1 + attr(num, "match.length"))
      perc <- as.numeric(perc)/100
      if (length(grep("^-.*$", calibration)) == 1) {
        ## calibration starts from distant (older) end
        cal_index <- 1:ceiling(perc * m)
        ver_index <- c(1:m)[-cal_index]
        cal_str <-
          gettextf("%d percent (= %d years) of data starting at older end",
                   perc * 100, length(cal_index))
      } else {
        ## calibration starts from recent end
        cal_index <- (m - floor(perc * m)):m
        ver_index <- c(1:m)[-cal_index]
        cal_str <-
          gettextf("%d percent (= %d years) of data starting at recent end",
                   perc * 100, length(cal_index))
      }
    } else {
      stop("Only percentage values can be supplied as character values to `calibration`, e.g. '66%' or '-45%'.")
    }
  } else {
    if (is.numeric(calibration)) {
      if (length(calibration) == 1) {
        if (calibration < m) {
          if (sign(calibration) == -1) {
            ## calibration starts x years from older end
            calibration <- abs(calibration)
            cal_index <- 1:calibration
            ver_index <- c(1:m)[-cal_index]
            cal_str <- gettextf("%d least recent observations", calibration)
          } else {
            ## calibration starts x years from younger end
            cal_index <- (m - calibration + 1):m
            ver_index <- c(1:m)[-cal_index]
            cal_str <- gettextf("%d most recent observations", calibration)
          }
        } else {
          stop(gettextf("The provided data has only %d years - consider adapting `calibration` accordingly.", m))
        }
      } else {
        if (!all(calibration %in% all_years)) {
          stop("Given years do not correspond to available years in the supplied data.")
        } else {
          ## calibration is given as range of years
          cal_index <- c(1:m)[all_years %in% calibration]
          ver_index <- c(1:m)[-cal_index]
          cal_str <- gettextf("the timespan from %d to %d",
                              calibration[1], max(calibration))
        }
      }
    }
  }  
  
  ## divide data into sets
  y    <- chrono
  full <- data.frame(
    y = y,
    x = x
    )
  cal  <- data.frame(
    y = y[cal_index],
    x = x[cal_index]
    )
  ver  <- data.frame(
    y = y[ver_index],
    x = x[ver_index]
    )
  
  ## models
  lm_cal <- lmodel2(x ~ y, data = cal,
                    range.x = "interval",
                    range.y = "interval",
                    nperm = 100)
  
  lm_full <- lmodel2(x ~ y, data = full,
                     range.x = "interval",
                     range.y = "interval",
                     nperm = 100)

  r_cal <- lm_cal$r
  r_full <- lm_full$r

  ## extract model specs
  if (model == "ols") {
    specs_cal <- subset(lm_cal$regression.results, Method == "OLS")
    specs_full <- subset(lm_full$regression.results, Method == "OLS")
  } else {
    specs_cal <- subset(lm_cal$regression.results, Method == "RMA")
    specs_full <- subset(lm_full$regression.results, Method == "RMA")
  }
  intercept_cal <- specs_cal[[2]]
  slope_cal <- specs_cal[[3]]
  p_cal <- specs_cal[[5]]
  intercept_full <- specs_full[[2]]
  slope_full <- specs_full[[3]]
  p_full <- specs_full[[5]]
  
  predict_cal <- intercept_cal + slope_cal * cal$y
  predict_ver <- intercept_cal + slope_cal * ver$y
  
  RE <- reduction_of_error(ver$x, predict_ver, cal$x)
  CE <- coefficient_of_efficiency(ver$x, predict_ver, ver$x)
  DW <- dwtest(cal$x ~ predict_cal)
  
  model_lm <- list(
    call       = mf,
    target     = X$names[1],
    r.cal      = r_cal,
    r.full     = r_full,
    coef.cal   = c(intercept = intercept_cal, slope = slope_cal),
    coef.full  = c(intercept = intercept_full, slope = slope_full),
    p.cal      = p_cal,
    p.full     = p_full,
    RE         = RE,
    CE         = CE,
    DW         = DW,
    cal.model  = lm_cal,
    full.model = lm_full,
    cal.str    = cal_str,
    pred.cal   = predict_cal,
    pred.ver   = predict_ver,
    cal        = cal,
    ver        = ver,
    full       = full,
    years      = all_years,
    cal.years  = all_years[cal_index],
    ver.years  = all_years[ver_index]
    )
  class(model_lm) <- c("tc_skills", "list")
  model_lm
}
