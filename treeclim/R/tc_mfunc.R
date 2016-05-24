##' moving response and correlation function
##' 
##' see doc for dcc for details.
##' @param chrono tree-ring chronology
##' @param climate data.frame with climate parameters
##' @param boot which bootstrapping method should be used? (or "none")
##' @param sb logical: draw statusbar or not?
##' @param start_last logical: start with last (oldest) window?
##' @param win_size numeric: size of the moving in years
##' @param win_offset numeric: size of offset between moving windows
##' in years
##' @param ci numeric: p-level for confidence interval (must be in
##' c(0.1, 0.05, 0.01)
##' @param method character: method to be used (one of "response" or
##' "correlation")
##' @param p probability for rgeom, that determines distribution of
##' sampling blocks for stationary bootstrap scheme
##' @keywords internal
tc_mfunc <- function(chrono, climate, boot, sb, start_last,
                     win_size, win_offset, ci, method) {

  vnames <- climate$names
  pretty_names <- climate$pretty_names
  ## number of windows
  years <- as.numeric(rownames(climate$aggregate))
  nyears <- length(years)
  win_num <- floor((nyears - win_size + 1) / win_offset)
  if (win_num < 2) {
    stop(paste("Less than 2 windows. Consider a timespan greater than ",
               nyears, " or a win_size smaller than ", win_size, ".",
               sep = ""))
  }
  win_years_string <- character(win_num)
  windows <- 1:win_num

  ## initialize result matrices
  result_matrix_coef <- result_matrix_ci_upper <-
    result_matrix_ci_lower <- result_matrix_significant <-
      matrix(NA, ncol = win_num, nrow = dim(climate$aggregate)[2])
  
  if (sb)                            # initialize status bar (if TRUE)
    mpb <- txtProgressBar(min = 1,  max = win_num, style = 3)

  for (k in 1:win_num) {
    if (start_last) {
      series_subset_index <- ((nyears - ((k-1) * win_offset)) -
                              (win_size - 1)):(nyears - ((k-1) * win_offset))
    } else {
      series_subset_index <- (1 + ((k-1) * win_offset)):(1 + ((k-1) *
                                                              win_offset) +
                                                         (win_size - 1))
    }
    
    climate_win <- climate$aggregate[series_subset_index,]
    ## recover the original list structure
    climate_win_list <- list(
      aggregate = climate_win,
      names = climate$names,
      param = climate$param,
      months = climate$month,
      pretty_names = climate$pretty_names
      )
    chrono_win <- chrono[series_subset_index]

    if (method == "response") {
      window <- tc_response(chrono_win, climate_win_list,
                            ci = ci, boot = boot)
    } else {
      window <- tc_correlation(chrono_win, climate_win_list,
                               ci = ci, boot = boot)
    }
    
    result_matrix_coef[,k] <- window$result$coef
    result_matrix_ci_upper[,k] <- window$result$ci_upper
    result_matrix_ci_lower[,k] <- window$result$ci_lower
    result_matrix_significant[,k] <- window$result$significant
    win_years_string[k] <- paste(years[series_subset_index][1],
                                 years[series_subset_index][win_size],
                                 sep = "-")
    
    if (sb)                             # update status bar (if TRUE)
      setTxtProgressBar(mpb, k)
    
  }

  ## reorder output if necessary
  if (start_last) {
    result_matrix_coef <- result_matrix_coef[,win_num:1]
    result_matrix_ci_upper <- result_matrix_ci_upper[,win_num:1]
    result_matrix_ci_lower <- result_matrix_ci_lower[,win_num:1]
    result_matrix_significant <- result_matrix_significant[,win_num:1]
    win_years_string <- win_years_string[win_num:1]
  }
  
  out <- list(
    result = list()
  )
  
  out$result$coef <- data.frame(result_matrix_coef)
  colnames(out$result$coef) <- win_years_string
  rownames(out$result$coef) <- vnames
  out$result$ci_upper <- data.frame(result_matrix_ci_upper)
  colnames(out$result$ci_upper) <- win_years_string
  rownames(out$result$ci_upper) <- vnames
  out$result$ci_lower <- data.frame(result_matrix_ci_lower)
  colnames(out$result$ci_lower) <- win_years_string
  rownames(out$result$ci_lower) <- vnames
  out$result$significant <- data.frame(result_matrix_significant)
  colnames(out$result$significant) <- win_years_string
  rownames(out$result$significant) <- vnames
  out$result$pretty_names <- pretty_names
  
  if (sb)                               # close status bar (if TRUE)
    close(mpb)
  class(out$result) <- c("tc_mcoef", "list")
  out
}
