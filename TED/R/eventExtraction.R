#' Extract events from time series 
#' 
#' This function returns the starting and ending points of events according to the noise test results from a time series.
#' 
#' @param tests test p values from the noist tests for the subsequences.
#' @param w sliding window size.
#' @param alpha the significance level. When the noise test p value of the subsequence is smaller than this significance level,
#' it is a potential event. Default is 0.05.
#' @return a list consisting:
#'  
#' \item{start}{a vector consisting of starting points of events.}
#' 
#' \item{end}{a vector consisting of ending points of events.}
#' 
#' \item{tests}{smoothed test \code{p} value series.}
#' 
#' \item{nevents}{number of detected events.}
#' @references Yanfei Kang, Danijel Belusic, Kate Smith-Miles (2014): Detecting and Classifying Events in Noisy Time Series. 
#' \emph{J. Atmos. Sci.}, \bold{71}, 1090-1104.
#' \url{http://dx.doi.org/10.1175/JAS-D-13-0182.1}.

#' @export

eventExtraction <- function(tests, w, alpha = 0.05) {
    
    N = length(tests)
    if (sum(is.na(tests)) == N) {
        warning("All values are missing for the current time series.")
        results <- list(start = c(), end = c(), tests = c(), nevents = 0)
        return(results)
    } else {
        values = rep(NA, N)
        values[which(tests > alpha & tests < 1)] = 2
        values[which(tests == 1)] = 1
        values[which(tests <= alpha)] = 0
        rleevents = rle(values)
        index1 = which((rleevents$values == 1))
        tests_new = tests
        # smoothing
        i = 1
        while (i <= length(index1) - 1) {
            if (sum(rleevents$length[(index1[i] + 1):(index1[i + 1] - 1)] <= ceiling(w/10)) == length(rleevents$length[(index1[i] + 
                1):(index1[i + 1] - 1)])) {
                tests_new[(sum(rleevents$lengths[1:(index1[i])]) + 1):(sum(rleevents$lengths[1:(index1[i + 1])]))] = 1
            }
            i = i + 1
        }
        
        values[which(tests_new > alpha & tests_new < 1)] = 2
        values[which(tests_new == 1)] = 1
        values[which(tests_new <= alpha)] = 0
        rleevents = rle(values)
        index0 = which((rleevents$values == 0) & (rleevents$lengths >= w/6))
        # smoothing
        i = 1
        while (i <= length(index0) - 1) {
            if (sum(rleevents$length[(index0[i] + 1):(index0[i + 1] - 1)] <= ceiling(w/10)) == length(rleevents$length[(index0[i] + 
                1):(index0[i + 1] - 1)])) {
                tests_new[(sum(rleevents$lengths[1:(index0[i])]) + 1):(sum(rleevents$lengths[1:(index0[i + 1])]))] = 0
            }
            if ((1 %in% rleevents$values[(index0[i] + 1):(index0[i + 1] - 1)][which(rleevents$length[(index0[i] + 1):(index0[i + 
                1] - 1)] > w/20)]) & (!1 %in% rleevents$values[(index0[i] + 1):(index0[i + 1] - 1)][which(rleevents$length[(index0[i] + 
                1):(index0[i + 1] - 1)] > w)]) & (sum(rleevents$length[(index0[i] + 1):(index0[i + 1] - 1)][which(rleevents$values[(index0[i] + 
                1):(index0[i + 1] - 1)] == 1)]) < w) & (!2 %in% rleevents$values[(index0[i] + 1):(index0[i + 1] - 1)][which(rleevents$length[(index0[i] + 
                1):(index0[i + 1] - 1)] > w/20)]) & (length(which(rleevents$length[(index0[i] + 1):(index0[i + 1] - 1)] < 5)) < 
                3)) {
                tests_new[(sum(rleevents$lengths[1:(index0[i])]) + 1):(sum(rleevents$lengths[1:(index0[i + 1])]))] = 0
            }
            i = i + 1
        }
        
        # defining potential events
        events = as.numeric(tests_new < alpha)
        # potential events should be sequentially enough to indicate the event
        rleevents = rle(events)
        lengths = rleevents$lengths
        eventslengths = lengths[(lengths > w/2) & (rleevents$values == 1)]
        
        index = which((lengths > w/2) & (rleevents$values == 1))
        if (length(index) == 0) {
            results <- list(start = c(), end = c(), tests = tests_new, nevents = 0)
            return(results)
        } else {
            points = rep(NA, length(index))
            for (i in 1:length(points)) {
                points[i] = sum(lengths[1:index[i]])
            }
            start = points - eventslengths + 1
            end = points
            results <- list(start = start, end = end, tests = tests_new, nevents = length(start))
            return(results)
        }
    }
    
    
} 
