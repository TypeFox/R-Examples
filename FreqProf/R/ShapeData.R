#' Convert data to moving sum/prop.
#' 
#' @param data.behavior a data.frame containing occurrence/nonoccurrence data in
#'   binary (0/1) format
#' @param window the window length to use in computing a moving sum or
#'   proportion
#' @param step the number of bins of which the data will be translated.
#' @param resolution the number of points contained in a bin
#' @param which giving the moving function to apply: sum or proportion
#' @return The data in a \code{freqprof} object.
#' @export
#' @examples
#' data(s58)
#' freqprof(s58)
freqprof = function(data.behavior,
                    window     = round(0.25 * nrow(data.behavior)),
                    step       = 1,
                    resolution = 1,
                    which      = c('sum','proportion')) {
  # selecting the appropriate moving function, according to the variable 'which'
  # by default, which is 'sum'
  if(length(which) > 1) which = 'sum'
  
  if(!(which %in% c('sum','proportion'))) {
    stop("possible values for variable 'which' are c('sum','proportion').")
  }
  
  # computing frequency profile
  freqprof = as.data.frame(apply(data.behavior, 2, 
                                 function(x) {
                                   movfun(x,
                                          n   = window,
                                          s   = step,
                                          r   = resolution,
                                          fun = which)$movfun
                                 })
                           )
  
  res = cbind(data.frame(time   = (0:(nrow(freqprof) - 1)) * resolution,
                         panels = movfun(data.behavior[, 1], 
                                         n   = window, 
                                         s   = step, 
                                         r   = resolution, 
                                         fun = which)$panels), freqprof)
  
  return(structure(list(window     = window,
                        step       = step,
                        resolution = resolution,
                        raw.data   = data.behavior,
                        type       = which,
                        data       = res),
                   class = "freqprof"))
}

#' Internal function for Resolution Adjustment
#' 
#' Internal function in \code{\link{freqprof}} 
#' that is used to modify data resolution.
#' 
#' @param x data data passed from \code{\link{freqprof}}
#' @param r resolution passed from \code{\link{freqprof}}
#' @return Resolution adjustment.
radj <- function(x, r) {
  # x is data
  # r is resolution
  adj <- rep(NA, floor((length(x) / r)))
  
  for(j in 1:length(adj)) {
    adj[j] <- sum(x[(1 + (j - 1) * r):(j * r)])
    adj[j] <- ifelse(adj[j] > 0, 1, 0)
  }
  
  return(adj)
}

#' Internal function for Generating Moving Sum or Proportion
#' 
#' Internal function in \code{\link{freqprof}} that is used to generate moving 
#' sum or proportion data.
#' 
#' @param x data passed from \code{\link{freqprof}}
#' @param n window length passed from \code{\link{freqprof}}
#' @param s step size passed from \code{\link{freqprof}}
#' @param r resolution passed from \code{\link{freqprof}}
#' @param fun "sum" or "proportion" passed from \code{\link{freqprof}}
#' @return Returns a list containing the processed data into $movfun, and the 
#'   associated panels into $panels. Passes list to \code{\link{freqprof}}.
movfun = function(x, n, s, r, fun) {
  if (r > 1) {
    x <- radj(x, r)
  }
  
  fun = switch(fun,
               sum        = sum,
               proportion = function(y) {sum(y) / n})
  
  res = rep(NA, floor((length(x) + n - 1) / s))
  
  panels = rep(2, floor((length(x) + n - 1) / s))
  
  for(j in 1:length(res)) {
    res[j] = fun(x[max(1, j * s - n + 1):min(j * s, length(x))])
    
    if((j * s - n + 1) < 1) panels[j] = 1
    if((j * s) > length(x)) panels[j] = 3
  }
  
  return(list(movfun = c(0, res, 0), 
              panels = c(1, panels, 3)))
}
