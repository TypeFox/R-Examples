kristrom <- function(formula, data, subset){
  
  if(missing(data)) data <- environment(formula)

  mf <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "subset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- formula
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  original.data <- data
  data <- mf

  # removing observations with missing values
  na.num <- max(sum(as.numeric(is.na(data))))
  if(na.num != 0){ 
    d1 <- nrow(data)
    data <- na.omit(data)
    d2 <- nrow(data)
    warning(paste("Missing values detected.", d1 - d2, "rows are removed.", sep = " "))
  }
  
  lhs <- formula[[2]]      # extracting the name of the acceptance/rejection variable from the formula supplied
  rhs <- formula[[3]]      # extracting the name of the bid variable from the formula supplied
  R <- eval(lhs, data)     # yes/no to the bids
  bid <- eval(rhs, data)   # the suggested bid

  nobs <- length(R)         # the number of observations
  data.table <- as.data.frame(cbind(bid, R)) # makig a data frame object with the bid and yes/no variables

  # making a matrix describing the number of respondents who answered 
  # "yes" to CV question, total number of respondents, and ratio of 
  # respondents who answered "yes" among the total number of respondents 
  # for each value of the suggested bids.
  tab.dat <- table(data.table)
  tab.dat <- cbind(tab.dat, rowSums(tab.dat))
  tab.dat <- cbind(tab.dat, tab.dat[,2]/tab.dat[,3])
  tab.dat <- tab.dat[, -1]
  colnames(tab.dat) <- c("yes", "total", "yes/total")  # adding column names
  tab.dat <- rbind(c(1, 1, 1), tab.dat)   # adding "observatin" for bid = 0
  rownames(tab.dat)[1] <- "0"
  unq.bid <- as.numeric(rownames(tab.dat))

  M <- nrow(tab.dat)  # the number of strata
  p <- tab.dat[,3]    # preliminary probability estimates
  adj.P <- p
  diff.p <- diff(adj.P)
  j <- which(diff.p > 0)[1]
  while(!is.na(j)){ # adjusting probabilities so that the survival function is non-increasing.
    tab.dat[j, 1] <- tab.dat[j+1, 1] <- tab.dat[j, 1] + tab.dat[j+1, 1]
    tab.dat[j, 2] <- tab.dat[j+1, 2] <- tab.dat[j, 2] + tab.dat[j+1, 2]
    tab.dat[, 3]  <- tab.dat[, 1]/tab.dat[, 2]
    adj.P <- tab.dat[, 3]
    diff.p <- diff(adj.P)
    j <- which(diff.p > 0)[1]
  }

#  estimates <- cbind(c(NA, unq.bid), c(unq.bid, Inf), c(adj.P, 0))
#  colnames(estimates) <- c("Lower", "Upper", "Prob.")
#  rownames(estimates) <- seq(1, nrow(estimates))

  # organizing a table with upper bounds and respective adjusted pribabilities
  estimates <- cbind(c(unq.bid, Inf), c(adj.P, 0))
  colnames(estimates) <- c("Upper", "Prob.")
  rownames(estimates) <- seq(1, nrow(estimates))

  # arranging outcomes into a single list variable
  output <- list(
        adj.P = adj.P, 
        M = M,        # the number of strata
        nobs = nobs, 
        unq.bid = unq.bid,
        tab.dat = tab.dat,
        estimates = estimates
      )

  class(output) <- "kristrom"
  return(output)

}

print.kristrom <- function(x, digits = 4, ...){

  cat("\nProbability:", "\n", sep = " ")
  print(x$estimates, digits = 4)
  invisible(x)

}

summary.kristrom <- function(object, digits = max(3, getOption("digits") - 1), ...){
  # the area under the empirical survival function up to the max bid
  area <- sum(0.5*(object$adj.P[-object$M] + object$adj.P[-1])*diff(object$unq.bid))
  # X-intercepr
  if(object$adj.P[object$M-1] != object$adj.P[object$M]){
    x.icpt <- object$unq.bid[object$M] + object$adj.P[object$M]*(object$unq.bid[object$M]-object$unq.bid[object$M-1])/
              (object$adj.P[object$M-1]-object$adj.P[object$M])
  } else if(object$adj.P[object$M-2] != object$adj.P[object$M-1]){
    x.icpt <- object$unq.bid[object$M] + object$adj.P[object$M]*(object$unq.bid[object$M-1]-object$unq.bid[object$M-2])/
              (object$adj.P[object$M-2]-object$adj.P[object$M-1])
  } else {
    x.icpt <- 1.5*object$unq.bid[object$M]    # this is somewhat arbitrary...
  }
  # mean WTP
  object$med.meanWTP <- area + 0.5*(x.icpt - object$unq.bid[object$M])*object$adj.P[object$M]         # treatment of the corner point?
  names(object$med.meanWTP) <- "meanWTP(SK)"

  object$meanWTP <- sum(object$adj.P[-1]*diff(object$unq.bid))         
  names(object$meanWTP) <- "meanWTP(KM)"       # the end point of the survival function calculated by linear interpolation
  object$x.icpt <- x.icpt
  names(object$x.icpt) <- "x intercept"
  # median estimation
  m <- max(which(object$adj.P > 0.5))
  ratio <- c(object$adj.P[m] - 0.5, 0.5 - object$adj.P[m+1])
  object$medianWTP <- object$unq.bid[m] + (object$unq.bid[m+1] - object$unq.bid[m])*ratio[1]/sum(ratio)
  names(object$medianWTP) <- "medianWTP"

  class(object) <- "summary.kristrom"
  return(object)

}

print.summary.kristrom <- function(x, digits = max(3, getOption("digits") - 1), ...){
  cat("Survival probability:", "\n", sep = " ")
  print.default(x$estimates, digits = 4, right = TRUE, print.gap = 2)

  cat("\nWTP estimates:", sep = " ")
  cat("\n Mean:", formatC(x$meanWTP, format="f", digits = digits), "", sep = " ")
  cat(" (Kaplan-Meier)", sep = "")
  cat("\n Mean:", formatC(x$med.meanWTP, format="f", digits = digits), "", sep = " ")
  cat(" (Spearman-Karber)", sep = "")
  # cat("\nMedian WTP in:", "[", formatC(x$medianWTP[1], digits = digits), ",", formatC(x$medianWTP[2], digits = digits), "]", "\n", sep = "")
  cat("\n Median:", formatC(x$medianWTP, format="f", digits = digits), "\n", sep = " ")

}

plot.kristrom<- function(x, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, lwd = NULL, lty = NULL, ...){

  plot.x <- summary.kristrom(x)  # summarizing the object for plot
  pr <- plot.x$estimates[, 2]    # probability estimates
  x.ax <- c(plot.x$unq.bid, plot.x$x.icpt)  # x-axis
  meanWTP <- plot.x$meanWTP
  medianWTP <- plot.x$medianWTP
  if(is.null(main)) main <- ""                      # main title
  if(is.null(sub)) sub <- ""                        # subtitle
  if(is.null(xlab)) xlab <- "Bid"                   # label of x-axis
  if(is.null(ylab)) ylab <- "Survival Probability"  # label of y-axis
  if(is.null(lwd)) lwd <- 3                         # line width
  if(is.null(lty)) lty <- 1                         # line type
  
  plot(x.ax, pr, axes = F, xlab = xlab, ylab = ylab, main = main, lwd = lwd, type = "l")
  axis(1, pos = 0, at = round(x.ax, 0), adj = 0)   # adding the x-axis
  axis(2, pos = 0, at = seq(0, 1, by = 0.2), las = 2, adj = 1)  # adding the y-axis
#  if(mean) segments(x0 = meanWTP, y0 = 0, y1 = pr[max(which(x.ax < meanWTP))], lty = 2)
#  if(median) segments(x0 = medianWTP, y0 = 0, y1 = pr[max(which(x.ax < medianWTP))], lty = 3)
#  if(median) segments(x0 = 0, x1 = x.ax[max(which(pr > 0.5))], y0 = 0.5, lty = 3)
}
