### Predictive Probability
glarmaPredProb <- function(object){
  if(object$type == "Poi"){
    ## Getting lambda from the model
    lambda <- fitted(object)
    ## Lower bound for conditional PIT: F(y-1)
    lower <- ppois(object$y - 1, lambda)
    ## Upper bound for conditional PIT: F(y)
    upper <- ppois(object$y, lambda)
  }
  if(object$type == "Bin"){
    pi <- exp(object$W) /
          (1 + exp(object$W))
    ## Number of trials
    nTrial <- rowSums(object$y)
    ## Lower bound for conditional PIT: F(y-1)
    lower <- pbinom(object$y[, 1] - 1, nTrial, pi)
    ## Upper bound for conditional PIT: F(y)
    upper <- pbinom(object$y[, 1] , nTrial, pi)
  }
  if(object$type == "NegBin"){
    ## Getting E(y), mu, from the model
    mu <- fitted(object)
    ## Gamma scale parameter estimated from the model
    alpha <- coef(object, "NB")
    ## E(y) = [alpha * (1 - prob)] / prob
    prob <- alpha / (mu + alpha)
    ## Lower bound for conditional PIT: F(y-1)
    lower <- pnbinom(object$y - 1, alpha, prob)
    ## Upper bound for conditional PIT: F(y)
    upper <- pnbinom(object$y, alpha, prob)
  }
  list("lower" = lower, "upper" = upper)
}

### Non-randomized PIT
glarmaPIT <- function(object, bins = 10) {
  ## Setting up the inputs
  dummy.variable <- seq(0, 1, by = 1 / bins)
  predProb <- glarmaPredProb(object)
  ## Remove the column of the failure in y for binomial distribution
  if (object$type == "Bin") object$y <- object$y[, 1]
  ## Computing the PIT conditioned on each y and its corresponding parameters
  con.PIT <- matrix(0, ncol = length(object$y), nrow = length(dummy.variable))
  for (i in 1:length(object$y)){
    for (j in 1:length(dummy.variable)){
      if(dummy.variable[j] <= predProb$lower[i]) con.PIT[j, i] = 0

      if(dummy.variable[j] > predProb$lower[i] &
           dummy.variable[j] < predProb$upper[i]){
        con.PIT[j, i] = (dummy.variable[j] - predProb$lower[i]) /
                        (predProb$upper[i] - predProb$lower[i])
      }
      if(dummy.variable[j] >= predProb$upper[i]) con.PIT[j, i] = 1
    }
  }
  ## Computing PIT. PIT = sum(Conditional PIT) / (T - 1)
  ## The column for PIT is the average of the column and
  ## the previous column for conditional PIT
  PIT <- matrix(0, ncol = length(object$y), nrow = length(dummy.variable))
  for (i in 2:length(object$y)){
    PIT[, i] <- apply(con.PIT[, 1:i], 1, sum) * (i - 1)^(-1)
  }
  list("ConditionalPIT" = con.PIT, "PIT" = PIT)
}

### PIT histogram
histPIT <- function(object, bins = 10, line = TRUE, colLine = "red",
                    colHist = "royal blue", lwdLine = 2,
                    main = NULL, ...) {
  ## Computing the height for the bar
  PIT <- glarmaPIT(object, bins = bins)$PIT
  height <- diff(PIT[, ncol(PIT)]) * bins
  ## Control ylim to show the peak values
  if (max(height) > 2){
    y.upper <- max(height) + 1 / (bins / 2)
  }
  else {
    y.upper <- 2
  }
  ## Plotting the bar graph
  if (is.null(main)) {
      main <- paste("PIT for GLARMA",
                    paste("(",
                          switch(object$type,
                                 "Poi" = "Poisson",
                                 "Bin" = "Binomial",
                                 "NegBin" = "Negative Binomial"),
                          ")", sep = ""))
  }
  barplot(height, ylim = c(0, y.upper),
          border = TRUE, space = 0, xpd = FALSE,
          xlab = "Probability Integral Transform",
          ylab = "Relative Frequency",
          main = main , col = colHist, ...)
  ## Plotting the comparison line
  if (line == TRUE){
    abline(h = 1, lty = 2, col = colLine, lwd = lwdLine)
  }
  ## Adjusting the x-axis scale
  plot.window(xlim = c(0, 1), ylim = c(0, y.upper))
  axis(1)
  box()
}

### Q-Q plot for PIT
qqPIT <- function(object, bins = 10, col1 = "red", col2 = "black",
                  lty1 = 1, lty2 = 2, type = "l",
                  main = NULL, ...){
  dummy.variable <- seq(0, 1, by = 1 / bins)
  PIT <- glarmaPIT(object, bins = bins)$PIT
  qq.plot <- PIT[, ncol(PIT)]
  if (is.null(main)){
      main <- paste("Uniform Q-Q plot",
                    paste("(", switch(object$type,
                                      "Poi" = "Poisson",
                                      "Bin" = "Binomial",
                                      "NegBin" = "Negative Binomial"),
                          ")", sep = ""))
  }
  plot(dummy.variable, qq.plot, lty = lty1, col = col1,
       xlim = c(0, 1), ylim = c(0, 1), type = type,
       xlab = "Theoretical", ylab = "Sample",
       main = main, ...)
  abline(0, 1, col = col2, lty = lty2)
  list("sample" = qq.plot, "theoretical" = dummy.variable)
  invisible()
}
