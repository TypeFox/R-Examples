#' Information value of an independent variable in predicting a binary response
#' 
#'  Takes in independent and dependent variable and returns IV value
#'  @param x an independent variable
#'  @param y a binary response variable
#'  @details
#'  Information value of a variable is a significant indicator of its relevance in
#'  the prediction of a binary response variable. \code{iv} computes that value using
#'  the formula,
#'  IV = summation[(Responders - Non-responders)*ln(Responders / Non-responders) for
#'  each bin].
#'  
#'  Ten bins are created for continous variables while categories itself are used as
#'  bins for categorical independent variables.
#'  @return information value of \code{x}
#'  @author Akash Jain
#'  @seealso \code{\link{accuracy}}, \code{\link{auc}}, \code{\link{ks}}, \code{\link{splitdata}}
#'  @examples
#'  # A 'data.frame'
#' df <- data.frame(x = c('a', 'a', 'a', 'b', 'b', 'b'),
#'                  y = c(0, 1, 0, 1, 0, 1))
#'
#' # Information value
#' IV <- iv(x = df[, 'x'], y = df[, 'y'])
#'  @export
iv <- function(x, y) {
  if(length(unique(y)) != 2 | (class(y) != 'integer' && class(y) != 'numeric' && class(y) != 'factor')) {
    stop('Invalid input: vector y should be integer or factor representing a binary response')
  } else if(length(x) != length(y)) {
    stop('Invalid input: vectors x and y should have the same length')
  } else if(class(x) != 'integer' && class(x) != 'numeric' && class(x) != 'character' && class(x) != 'factor') {
    stop('Invalid input: vector x should have class integer or numeric or character or factor')
  } else {
    if(class(x) == 'integer' | class(x) == 'numeric') {
      ntile <- as.integer(as.character(cut(x, 
                                           breaks = quantile(x, probs = seq(0, 1, by = 0.1), na.rm = TRUE), 
                                           labels = 1:10, 
                                           include.lowest = TRUE)))
    } else if(class(x) == 'character' | class(x) == 'factor') {
      ntile <- x
    }
    data <- data.frame(dec = ntile, y = y)
    dataNomiss <- data[!is.na(data$dec), ]
    dataNomiss$y <- as.integer(as.character(dataNomiss$y))
    t1 <- aggregate(y ~ dec, dataNomiss, function(x) sum(x == 0))
    t2 <- aggregate(y ~ dec, dataNomiss, function(x) sum(x == 1))
    summary <- merge(t1, t2, by = c('dec'), all.x = TRUE)
    names(summary)[2:3] <- c('numZero', 'numOne')
    sumZeroOv <- sum(as.character(y) == '0')
    sumOneOv <- sum(as.character(y) == '1')
    summary[, 'b'] <- summary[, 'numZero']/sumZeroOv
    summary[, 'g'] <- summary[, 'numOne']/sumOneOv
    summary[, 'logPart'] <- log(summary[, 'g']/summary[, 'b'])
    summary[ , 'IV'] <- (summary[, 'g'] - summary[, 'b']) * summary[, 'logPart']
    iv <- sum(summary[, 'IV'], na.rm = TRUE)
    return(iv)
  }
}