#' Kolmogorov-Smirnov statistic for predicted binary response
#' 
#'  Takes in actual binary response and predicted probabilities, and returns table 
#'  used for KS computation and KS value
#'  @param y actual binary response
#'  @param yhat predicted probabilities corresponding to the actual binary response
#'  @details
#'  Kolmogorov-Smirnov statistic can be easily computed using \code{ks} function. It
#'  not only computes the statistic but also returns the table used for ks computation.
#'  A chart is also plotted for quick visualization of split between percentage of
#'  responders and non-responders. Deciles are used for computation.
#'  
#'  Appropriate elements of the 2 element list (ie. ksTable or ks) can be picked
#'  using the code given in the example.  
#'  @return a two element list: table for KS computation and KS value itself
#'  @author Akash Jain
#'  @seealso \code{\link{accuracy}}, \code{\link{auc}}, \code{\link{iv}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame' with y and yhat
#' df <- data.frame(y = c(1, 0, 1, 1, 0),
#'                  yhat = c(0.86, 0.23, 0.65, 0.92, 0.37))
#'
#' # KS table and value
#' ltKs <- ks(y = df[, 'y'], yhat = df[, 'yhat'])
#' ksTable <- ltKs$ksTable
#' KS <- ltKs$ks
#'  @export
ks <- function(y, yhat) {
  if(length(unique(y)) != 2 | (class(y) != 'integer' && class(y) != 'numeric' && class(y) != 'factor')) {
    stop('Invalid input: y should be integer or factor vector representing a binary response')
  } else if(class(yhat) != 'numeric' | max(yhat) > 1 | min(yhat) < 0) {
    stop('Invalid input: yhat should be numeric vector of predicted probabilities in the range 0 to 1')
  } else if(length(y) != length(yhat)) {
    stop('Invalid input: vectors y and yhat should have the same length')
  } else {
    data <- data.frame(y = as.integer(as.character(y)), yhat = yhat)
    data[, 'dec'] <- decile(data[, 'yhat'], decreasing = TRUE)
    dataNomiss <- data[!is.na(data$dec), ]
    t1 <- aggregate(y ~ dec, dataNomiss, FUN = function(x) length(!is.na(x)))
    t2 <- aggregate(y ~ dec, dataNomiss, FUN = sum, na.rm = TRUE)
    ksTable <- merge(t1, t2, by = c('dec'), all.x = TRUE)
    names(ksTable)[2:3] <- c('total', 'responders')
    ksTable[, 'nonresponders'] <- ksTable[, 'total'] - ksTable[, 'responders']
    ksTable[, 'cumResponders'] <- cumsum(ksTable[, 'responders'])
    ksTable[, 'cumNonresponders'] <- cumsum(ksTable[, 'nonresponders'])
    ksTable[, 'perCumResponders'] <- round(ksTable[, 'cumResponders']/max(ksTable[, 'cumResponders']), digits = 2) * 100
    ksTable[, 'perCumNonresponders'] <- round(ksTable[, 'cumNonresponders']/max(ksTable[, 'cumNonresponders']), digits = 2) * 100
    ksTable[, 'split'] <- ksTable[, 'perCumResponders'] - ksTable[, 'perCumNonresponders']
    ks <- max(ksTable[, 'split'])
    lt <- list(ksTable = ksTable,
               ks = ks)
    par(xaxs = 'i', yaxs = 'i')
    plot(x = c(0, ksTable$dec), 
         y = c(0, ksTable$perCumResponders), 
         type = 'l', 
         xlab = 'Decile', 
         ylab = 'Cumulative %',
         xlim = c(0, 10),
         ylim = c(0, 100),
         col = 'blue',
         lwd = 2)
    lines(x = c(0, ksTable$dec), 
          y = c(0, ksTable$perCumNonresponders), 
          type = 'l',
          xlim = c(0, 10),
          ylim = c(0, 100),
          col = 'red',
          lwd = 2)
    return(lt)
  }
}