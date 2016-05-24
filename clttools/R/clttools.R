#' Theoretical Probability Distribution of General Experiment
#'
#' General experiment with basic probability
#' @usage expt(x, n, prob = NULL)
#' @param x vector, possible outcomes in one trial of experiment
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples expt(x = c(1:3), n = 4)
#' expt(c(2:4), 3, prob = c(0.3, 0.5, 0.2))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

expt <- function (x, n, prob = NULL)
{
  sample_space <- expand.grid(rep(list(x),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                length(sample_space[,1]))
  }
  else {
    if (!identical(length(x), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }

    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }
  MEAN_VALUE <- apply(sample_space[,-(n+1)], 1, mean)
  PROBABILITY <- sample_space[, n+1]
  DF <- cbind(MEAN_VALUE, PROBABILITY)
  clt.mean <- aggregate(PROBABILITY~., data=DF, FUN = sum)
  return(clt.mean)
}

#' Theoretical Probability Distribution Plot of General Experiment
#'
#' General experiment plot with basic probability
#' @usage expt.plot(x, n, prob = NULL, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param x vector, possible outcomes in one trial of experiment
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples expt.plot(x = c(1:3), n = 4, col ='red', type = 'p')
#' expt.plot(c(2:4), 3, prob = c(0.3, 0.5, 0.2))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'


expt.plot <- function (x, n, prob = NULL, col = 'black',
                             type = NULL, main = NULL , sub = NULL)
{
  sample_space <- expand.grid(rep(list(x),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                length(sample_space[,1]))
  }
  else {
    if (!identical(length(x), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }

    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist, KEEP.OUT.ATTRS = FALSE)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }
  clt.mean <- tapply(sample_space[, n+1],
                      apply(sample_space[,-(n+1)], 1, mean),
                        FUN = sum)
  x <- as.numeric(names(clt.mean))
  y <- as.numeric(clt.mean)
  par(mfrow=c(1,1))
  plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = main , sub = sub)
}

#' Theoretical Probability Distribution of Rolling Dice
#'
#' Mean and probability of rolling fair or loaded dice
#' @usage dice(n, prob = NULL)
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples dice(n = 4)
#' dice(2, c(0.1, 0.2, 0.2, 0.1, 0.3, 0.1))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

dice <- function(n, prob = NULL){
  sample_space <- expand.grid(rep(list(c(1:6)),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                length(sample_space[,1]))
  }
  else {
    if (!identical(length(c(1:6)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }
  MEAN_VALUE <- apply(sample_space[,-(n+1)], 1, mean)
  PROBABILITY <- sample_space[, n+1]
  DF <- cbind(MEAN_VALUE, PROBABILITY)
  clt.mean <- aggregate(PROBABILITY~., data=DF, FUN = sum)
  return(clt.mean)
}

#' Theoretical Probability Distribution Plot of Rolling Dice
#'
#' Probability plot of rolling fair or loaded dice
#' @usage dice.plot(n, prob = NULL, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples dice.plot(n = 4, col ='red', type = 'p')
#' dice.plot(3, prob = c(0.3, 0.1, 0.2, 0.1, 0.1, 0.2))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

dice.plot <- function (n, prob = NULL, col = 'black', type = NULL,
                       main = NULL , sub = NULL)
{
  sample_space <- expand.grid(rep(list(c(1:6)),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                 length(sample_space[,1]))
  }
  else {
    if (!identical(length(c(1:6)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist, KEEP.OUT.ATTRS = FALSE)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }
  clt.mean <- tapply(sample_space[, n+1],
                     apply(sample_space[,-(n+1)], 1, mean), FUN = sum)
  x <- as.numeric(names(clt.mean))
  y <- as.numeric(clt.mean)
  par(mfrow=c(1,1))
  plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = main , sub = sub)
}


#' Theoretical Probability Distribution of Flipping Coins
#'
#' Mean and probability of flipping fair or loaded coin
#' @usage coin(n, prob = NULL)
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples coin(n = 4)
#' coin(6, c(0.1, 0.9))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

coin <- function(n, prob = NULL){
  sample_space <- expand.grid(rep(list(c(0:1)),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                 length(sample_space[,1]))
  }
  else {
    if (!identical(length(c(0:1)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }
  MEAN_VALUE <- apply(sample_space[,-(n+1)], 1, mean)
  PROBABILITY <- sample_space[, n+1]
  DF <- cbind(MEAN_VALUE, PROBABILITY)
  clt.mean <- aggregate(PROBABILITY~., data=DF, FUN = sum)
  return(clt.mean)
}

#' Theoretical Probability Distribution Plot of Flipping Coins
#'
#' Probability plot of flipping fair or loaded coin
#' @usage coin.plot(n, prob = NULL, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param n number of trials
#' @param prob probability assigned to each possible outcome
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all possible outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples coin.plot(n = 4, col ='red', type = 'p')
#' coin.plot(3, prob = c(0.3, 0.7))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

coin.plot <- function (n, prob = NULL, col = 'black', type = NULL,
                       main = NULL , sub = NULL)
{
  sample_space <- expand.grid(rep(list(c(0:1)),n))
  if (is.null(prob)) {
    sample_space$prob <- rep(1/length(sample_space[,1]),
                                 length(sample_space[,1]))
  }
  else {
    if (!identical(length(c(0:1)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
    probs <- expand.grid(plist, KEEP.OUT.ATTRS = FALSE)
    psample_space <- apply(probs, 1, prod)
    sample_space$prob <- psample_space
  }

  clt.mean <- tapply(sample_space[, n+1],
                     apply(sample_space[,-(n+1)], 1, mean), FUN = sum)
  x <- as.numeric(names(clt.mean))
  y <- as.numeric(clt.mean)
  par(mfrow=c(1,1))
  plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = main , sub = sub)
}


#' Probability Distribution of Simulated General Experiments
#'
#' Mean and probability of general simulated experiments
#' @usage expt.simu(x, n, times, prob = NULL)
#' @param x vector, possible outcomes in one trial of experiment
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples expt.simu(x = c(1:3), n = 4, times = 1000)
#' expt.simu(c(1:3), 4, 1000, prob = c(0.3, 0.1, 0.6))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

expt.simu <- function(x, n, times, prob= NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/length(x),length(x))
  }
  else {
    if (!identical(length(x), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  RM <- rowMeans(replicate(n, sample(x, size=times, replace=T, prob=prob)))
  MEAN_VALUE <- as.numeric(names(table(RM)))
  FREQUENCY <- table(RM)
  PROBABILITY <- FREQUENCY/times
  return(cbind(MEAN_VALUE,FREQUENCY,PROBABILITY))
}


#' Probability Distribution Plot of Simulated General Experiments
#'
#' Probability plot of general simulated experiments
#' @usage expt.simu.plot(x, n, times, prob = NULL, qqplot = FALSE, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param x vector, possible outcomes in one trial of experiment
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples expt.simu.plot(x = c(1:3), n = 4, times = 1000, col = 'red')
#' expt.simu.plot(c(1:3), 4, 1000, prob = c(0.3, 0.1, 0.6), type = 'p')
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

expt.simu.plot <- function (x, n, times, prob = NULL, qqplot = FALSE, col = 'black',
                                    type = NULL, main = NULL , sub = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/length(x),length(x))
  }
  else {
    if (!identical(length(x), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
   plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  if (qqplot == FALSE){
    par(mfrow=c(1,1))
    RM <- rowMeans(replicate(n, sample(x, size=times, replace=T, prob=prob)))
    MEAN_VALUE <- as.numeric(names(table(RM)))
    FREQUENCY <- table(RM)
    PROBABILITY <- FREQUENCY/times
    x <- MEAN_VALUE
    y <- PROBABILITY
    plot(x, y, xlab = "mean", ylab = "probability",  col = col,
         type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM <- rowMeans(replicate(n, sample(x, size=times, replace=T, prob=prob)))
    MEAN_VALUE <- as.numeric(names(table(RM)))
    FREQUENCY <- table(RM)
    PROBABILITY <- FREQUENCY/times
    x <- MEAN_VALUE
    y <- PROBABILITY
    plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
    qqnorm(RM)
    qqline(RM)}
}


#' Probability Distribution of Simulated Dice Rolling
#'
#' Mean and probabilityf of flipping fair or loaded dice
#' @usage dice.simu(n, times, prob = NULL)
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples dice.simu(n = 4, times = 1000)
#' dice.simu(4, 1000, prob = c(0.3, 0.1, 0.1, 0.1, 0.3, 0.1))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

dice.simu <- function(n, times, prob = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/6,6)
  }
  else {
    if (!identical(length(c(1:6)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist = list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }

  RM <- rowMeans(replicate(n, sample(c(1:6), size=times,
                                    replace=T, prob=prob)))
  MEAN_VALUE <- as.numeric(names(table(RM)))
  FREQUENCY <- table(RM)
  PROBABILITY <- FREQUENCY/times
  return(cbind(MEAN_VALUE,FREQUENCY,PROBABILITY))
}

#' Probability Distribution Plot of Simulated Dice Rolling
#'
#' Probability plot of dice simulated experiments
#' @usage dice.simu.plot(n, times, prob = NULL, qqplot = FALSE, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples dice.simu.plot(n = 4, times = 1000, col = 'red')
#' dice.simu.plot(4, 1000, prob = c(0.3, 0.1, 0.1, 0.1, 0.1, 0.3), type = 'p')
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

dice.simu.plot <- function (n, times, prob = NULL, qqplot = FALSE, col = 'black', type = NULL,
                            main = NULL , sub = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/6,6)
  }
  else {
    if (!identical(length(c(1:6)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM <- rowMeans(replicate(n, sample(c(1:6), size=times,
                                    replace=T, prob=prob)))
  MEAN_VALUE <- as.numeric(names(table(RM)))
  FREQUENCY <- table(RM)
  PROBABILITY <- FREQUENCY/times
  x <- MEAN_VALUE
  y <- PROBABILITY
  plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
  }

  if (qqplot == TRUE) {
    par(mfrow=c(1,2))
    RM <- rowMeans(replicate(n, sample(c(1:6), size=times,
                                       replace=T, prob=prob)))
    MEAN_VALUE <- as.numeric(names(table(RM)))
    FREQUENCY <- table(RM)
    PROBABILITY <- FREQUENCY/times
    x <- MEAN_VALUE
    y <- PROBABILITY
    plot(x, y, xlab = "mean", ylab = "probability",  col = col,
         type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Probability Distribution of Simulated Coins Flipping
#'
#' Mean and probability plot of flipping fair or loaded coin
#' @usage coin.simu(n, times, prob = NULL)
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @return Mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples coin.simu(n = 4, times = 1000)
#' coin.simu(4, 1000, prob = c(0.3, 0.7))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

coin.simu <- function(n, times, prob = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/2,2)
  }
  else {
    if (!identical(length(c(0:1)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  RM <- rowMeans(replicate(n, sample(c(0:1), size=times,
                                    replace=T, prob=prob)))
  MEAN_VALUE <- as.numeric(names(table(RM)))
  FREQUENCY <- table(RM)
  PROBABILITY <- FREQUENCY/times
  return(cbind(MEAN_VALUE,FREQUENCY,PROBABILITY))
}


#' Probability Distribution Plot of Simulated Coins Flipping
#'
#' Probability plot of simulated experiments on flipping coins
#' @usage coin.simu.plot(n, times, prob = NULL, qqplot = FALSE, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples coin.simu.plot(n = 4, times = 1000, col = 'red')
#' coin.simu.plot(4, 1000, prob = c(0.3, 0.7), type = 'p')
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

coin.simu.plot <- function (n, times, prob = NULL, qqplot = FALSE, col = 'black', type = NULL,
                            main = NULL, sub = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/2,2)
  }
  else {
    if (!identical(length(c(0,1)), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM <- rowMeans(replicate(n, sample(c(0:1), size=times,
                                    replace=T, prob=prob)))
  MEAN_VALUE <- as.numeric(names(table(RM)))
  FREQUENCY <- table(RM)
  PROBABILITY <- FREQUENCY/times
  x <- MEAN_VALUE
  y <- PROBABILITY
  plot(x, y, xlab = "mean", ylab = "probability",  col = col,
       type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
  }

  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM <- rowMeans(replicate(n, sample(c(0:1), size=times,
                                       replace=T, prob=prob)))
    MEAN_VALUE <- as.numeric(names(table(RM)))
    FREQUENCY <- table(RM)
    PROBABILITY <- FREQUENCY/times
    x <- MEAN_VALUE
    y <- PROBABILITY
    plot(x, y, xlab = "mean", ylab = "probability",  col = col,
         type = type, main = bquote(paste("n=",.(n), ', times=', .(times))) , sub = sub)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Poisson distribution
#'
#' Histogram and Q-Q plot of simulated Poisson distribution
#' @usage pois.simu.plot(n, lambda, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param lambda parameter of Poisson distribution
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Poisson distribution, red curve represents theoretical density
#' @examples pois.simu.plot(n = 5, lambda = 3, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

pois.simu.plot <- function(n, lambda, times,  ylim = NULL, qqplot = FALSE)
{
  lambda = lambda
    if (qqplot == FALSE){
    par(mfrow = c(1,1))
    RM = apply(matrix(rpois(n*times,lambda),nrow=times), 1, mean)
    hist(RM, main = "Histogram",freq=FALSE, xlab = "mean",  ylim = ylim)
    x = seq(min(RM),max(RM),length = 100)
    curve(dnorm(x,lambda,sd = sqrt(lambda/n)),col = "red", add = TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow = c(1,2))
    RM = apply(matrix(rpois(n*times,lambda),nrow=times), 1, mean)
    hist(RM, main = "Histogram",freq=FALSE, xlab = "mean",  ylim = ylim)
    x = seq(min(RM),max(RM),length = 100)
    curve(dnorm(x,lambda,sd = sqrt(lambda/n)),col = "red", add = TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Binomial distribution
#'
#' Histogram and Q-Q plot of simulated Binomial distribution
#' @usage binom.simu.plot(n, size, prob, times, ylim = NULL, qqplot = FALSE)
#' @param n number of observations
#' @param size number of trials (zero or more)
#' @param prob probability of success on each trial
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Binomial distribution, red curve represents theoretical density
#' @examples binom.simu.plot(n = 10, size = 5, prob = 0.2, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

binom.simu.plot <- function(n, size, prob, times, ylim = NULL, qqplot = FALSE)
{
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM = apply(matrix(rbinom(n*times, size = size, prob = prob),nrow=times), 1, mean)
  hist(RM, main = "Histogram", freq=FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x,size*prob,sd=sqrt(size*prob*(1-prob)/n)),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM = apply(matrix(rbinom(n*times, size = size, prob = prob),nrow=times), 1, mean)
    hist(RM, main = "Histogram", freq=FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x,size*prob,sd=sqrt(size*prob*(1-prob)/n)),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Negative Binomial distribution
#'
#' Histogram and Q-Q plot of simulated Negative Binomial distribution
#' @usage nbinom.simu.plot(n, size, prob, times, ylim = NULL, qqplot = FALSE)
#' @param n number of observations
#' @param size number of trials (zero or more)
#' @param prob probability of success on each trial
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Negative Binomial distribution, red curve represents theoretical density
#' @examples nbinom.simu.plot(n = 10, size = 5, prob = 0.2, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

nbinom.simu.plot <- function(n, size, prob, times, ylim = NULL, qqplot = FALSE)
{
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM = apply(matrix(rnbinom(n*times, size = size, prob = prob),nrow=times), 1, mean)
  hist(RM, main = "Histogram",freq= FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, size*(1-prob)/prob, sd=sqrt(size*(1-prob)/(n*prob**2))),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM = apply(matrix(rnbinom(n*times, size = size, prob = prob),nrow=times), 1, mean)
    hist(RM, main = "Histogram",freq= FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, size*(1-prob)/prob, sd=sqrt(size*(1-prob)/(n*prob**2))),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Geometric distribution
#'
#' Histogram and Q-Q plot of simulated Geometric distribution
#' @usage geom.simu.plot(n, prob, times, ylim = NULL, qqplot = FALSE)
#' @param n number of observations
#' @param prob probability of success on each trial
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Geometric distribution, red curve represents theoretical density
#' @examples geom.simu.plot(n = 10, prob = 0.2, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

geom.simu.plot <- function(n, prob, times, ylim = NULL, qqplot = FALSE)
{
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM = apply(matrix(rgeom(n*times, prob = prob),nrow=times), 1, mean)
  hist(RM, main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, (1-prob)/prob, sd=sqrt((1-prob)/(n*prob**2))),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM = apply(matrix(rgeom(n*times, prob = prob),nrow=times), 1, mean)
    hist(RM, main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, (1-prob)/prob, sd=sqrt((1-prob)/(n*prob**2))),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Hypergeometric distribution
#'
#' Histogram and Q-Q plot of simulated Hypergeometric distribution
#' @usage hyper.simu.plot(n, a, b, k, times, ylim = NULL, qqplot = FALSE)
#' @param n number of observations
#' @param a the number of white balls in the urn
#' @param b the number of black balls in the urn
#' @param k the number of balls drawn from the urn
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Hypergeometric distribution, red curve represents theoretical density
#' @examples hyper.simu.plot(n = 10, a = 10, b = 10, k = 5, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

hyper.simu.plot <- function(n, a, b, k, times, ylim = NULL, qqplot = FALSE)
{
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  RM = apply(matrix(rhyper(n*times, m = a, n = b, k = k),nrow=times), 1, mean)
  hist(RM, main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  prob = a/(a+b)
  curve(dnorm(x, k*prob, sd=sqrt(k*prob*(1-prob)*(a+b-k)/(a+b-1)/n)),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    RM = apply(matrix(rhyper(n*times, m = a, n = b, k = k),nrow=times), 1, mean)
    hist(RM, main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    prob = a/(a+b)
    curve(dnorm(x, k*prob, sd=sqrt(k*prob*(1-prob)*(a+b-k)/(a+b-1)/n)),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}



#' Histogram and Q-Q plot of simulated Uniform distribution
#'
#' Histogram and Q-Q plot of simulated Uniform distribution
#' @usage unif.simu.plot(n, min = 0, max = 1, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param min possible minimum value of Uniform distribution. Must be finite
#' @param max possible maximum value of Uniform distribution. Must be finite
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Uniform distribution, red curve represents theoretical density
#' @examples unif.simu.plot(n = 5, min = 3, max = 5, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

unif.simu.plot <- function(n, min = 0, max = 1, times, ylim = NULL, qqplot = FALSE)
{
   RM = apply(matrix(runif(n*times, min = min, max = max), nrow=times), 1, mean)
   if (qqplot == FALSE){
   par(mfrow=c(1,1))
   hist(RM,  main = "Histogram", freq= FALSE, xlab = "mean", ylim = ylim)
   x=seq(min(RM),max(RM),length=100)
   curve(dnorm(x,(max + min)/2, sd=sqrt((max - min)**2/(12*n))),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
   par(mfrow=c(1,2))
   hist(RM,  main = "Histogram", freq= FALSE, xlab = "mean", ylim = ylim)
   x=seq(min(RM),max(RM),length=100)
   curve(dnorm(x,(max + min)/2, sd=sqrt((max - min)**2/(12*n))),col="red", add=TRUE)
   qqnorm(RM)
   qqline(RM)
  }
}

#' Histogram and Q-Q plot of simulated Exponential distribution
#'
#' Histogram and Q-Q plot of simulated Exponential distribution
#' @usage expo.simu.plot(n, rate = 1, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param rate vector of rates
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Exponential distribution, red curve represents theoretical density
#' @examples expo.simu.plot(n = 5, rate = 2, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

expo.simu.plot <- function(n, rate = 1, times, ylim = NULL, qqplot = FALSE )
{

  RM = apply(matrix(rexp(n*times, rate = rate), nrow=times), 1, mean)
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, 1/rate, sd=sqrt(1/(n*rate**2))),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, 1/rate, sd=sqrt(1/(n*rate**2))),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Chi-Squared distribution
#'
#' Histogram and Q-Q plot of simulated Chi-Squared distribution
#' @usage chisq.simu.plot(n, df, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param df degrees of freedom (non-negative, but can be non-integer)
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Chi-Squared distribution, red curve represents theoretical density
#' @examples chisq.simu.plot(n = 5, df = 4, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

chisq.simu.plot <- function(n, df, times, ylim = NULL, qqplot = FALSE )
{
  RM = apply(matrix(rchisq(n*times, df = df), nrow=times), 1, mean)
  if (qqplot == FALSE){
    par(mfrow=c(1,1))
  hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, df, sd=sqrt(2*df/n)), col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, df, sd=sqrt(2*df/n)), col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Gamma distribution
#'
#' Histogram and Q-Q plot of simulated Gamma distribution
#' @usage gamm.simu.plot(n, shape, rate = 1, scale = 1/rate, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param shape shape parameter
#' @param rate vector of rates
#' @param scale scale parameter
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Gamma distribution, red curve represents theoretical density
#' @examples gamm.simu.plot(n = 5, shape = 3, rate = 1, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

 gamm.simu.plot <- function(n, shape, rate = 1, scale=1/rate, times, ylim = NULL, qqplot = FALSE)
{
  RM = apply(matrix(rgamma(n*times, shape, rate = rate),
                    nrow=times), 1, mean)
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, shape/rate, sd=sqrt(shape/(n*rate**2))),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, shape/rate, sd=sqrt(shape/(n*rate**2))),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Beta distribution
#'
#' Histogram and Q-Q plot of simulated Beta distribution
#' @usage beta.simu.plot(n, shape1, shape2, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param shape1 non-negative parameters of the Beta distribution
#' @param shape2 non-negative parameters of the Beta distribution
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Beta distribution, red curve represents theoretical density
#' @examples beta.simu.plot(n = 5, shape1 = 3, shape2 = 1, times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

beta.simu.plot <- function(n, shape1, shape2, times, ylim = NULL, qqplot = FALSE)
{
  RM = apply(matrix(rbeta(n*times, shape1 = shape1, shape2 = shape2),
                    nrow=times), 1, mean)
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, shape1/(shape1+shape2),
              sd=sqrt(shape1*shape2/((shape1+shape2)**2*n*(shape1+shape2+1)))),
              col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, shape1/(shape1+shape2),
                sd=sqrt(shape1*shape2/((shape1+shape2)**2*n*(shape1+shape2+1)))),
          col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Histogram and Q-Q plot of simulated Normal distribution
#'
#' Histogram and Q-Q plot of simulated Normal distribution
#' @usage normal.simu.plot(n, mean=0, sd=1, times, ylim = NULL, qqplot = FALSE)
#' @param n number of trials in one simulation
#' @param mean vector of means
#' @param sd vector of standard deviations
#' @param times number of simulations
#' @param ylim range of y-axis
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @return Histogram and Q-Q plot of simulated Normal distribution, red curve represents theoretical density
#' @examples normal.simu.plot(n = 5, mean = 3, sd =2,  times = 100)
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

normal.simu.plot <- function(n, mean=0, sd=1, times, ylim = NULL, qqplot = FALSE)
{
  RM = apply(matrix(rnorm(n*times, mean = mean, sd = sd), nrow=times), 1, mean)
  if (qqplot == FALSE){
  par(mfrow=c(1,1))
  hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
  x=seq(min(RM),max(RM),length=100)
  curve(dnorm(x, mean, sd = sd/sqrt(n)),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean", ylim = ylim)
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, mean, sd = sd/sqrt(n)),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}


#' Mean square error of simulated experiments
#'
#' Mean square error of simulated experiments
#' @usage expt.mse(x, n, times, prob = NULL)
#' @param x vector, possible outcomes in one trial of experiment
#' @param n number of trials
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @return Mean square error of simulated experiments
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples expt.mse(x = c(1:3), n = 4, times = 100)
#' expt.mse(c(0.1, 4, 2), 3, times = 50, prob = c(0.3, 0.5, 0.2))
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif
#'

expt.mse <- function (x, n, times, prob = NULL )
{
  expt <- function (x, n, prob = prob)
  {
    sample_space <- expand.grid(rep(list(x),n))
    if (is.null(prob)) {
      sample_space$prob <- rep(1/length(sample_space[,1]), length(sample_space[,1]))
    }
    else {
      if (!identical(length(x), length(prob))) {
        stop("'prob' is not the same length as 'possible outcomes'")
      }
      if (any(prob < 0)) {
        stop("'prob' contains negative values")
      }
      plist <- list()
      for (i in 1:n) {
        plist[[i]] = prob
      }
      probs <- expand.grid(plist)
      psample_space <- apply(probs, 1, prod)
      sample_space$prob <- psample_space
    }
    MEAN_VALUE <- apply(sample_space[,-(n+1)], 1, mean)
    PROBABILITY <- sample_space[, n+1]
    DF <- cbind(MEAN_VALUE, PROBABILITY)
    clt.mean <- aggregate(PROBABILITY~., data=DF, FUN = sum)
    #return(clt.mean)
    prob.theory = clt.mean
  }
  expt.simu <- function(x, n, times, prob= prob)
  {
    if (is.null(prob)) {
      prob <- rep(1/length(x),length(x))
    }
    else {
      if (!identical(length(x), length(prob))) {
        stop("'prob' is not the same length as 'possible outcomes'")
      }
      if (any(prob < 0)) {
        stop("'prob' contains negative values")
      }
      plist <- list()
      for (i in 1:n) {
        plist[[i]] = prob
      }
    }
    RM <- rowMeans(replicate(n, sample(x, size=times, replace=T, prob=prob)))
    MEAN_VALUE <- as.numeric(names(table(RM)))
    FREQUENCY <- table(RM)
    PROBABILITY <- FREQUENCY/times
    prob.simu = cbind(MEAN_VALUE,PROBABILITY)
  }
  prob.theory = expt(x, n, prob = prob)
  prob.simu = expt.simu(x, n, times, prob= prob)
  DF= merge(prob.theory, prob.simu, by="MEAN_VALUE", all= TRUE)
  DF[is.na(DF)] <- 0
  MSE = mean((DF[,2]-DF[,3])**2)
  return(MSE)
}



#' Histogram and Q-Q plot of any given continuous distribution
#'
#' Histogram and Q-Q plot of any given continuous distribution
#' @usage distr.simu.plot(distr, n, times, prob = NULL, qqplot = FALSE, col = "black", type = NULL,
#' main = NULL, sub = NULL)
#' @param distr vector, all possible outcomes, population distribution
#' @param n number of trials in one simulation
#' @param times number of simulations
#' @param prob probability assigned to each possible outcome
#' @param qqplot an argument to output Q-Q plot or not, can be TRUE or FALSE
#' @param col color of the plot
#' @param type type of plot
#' @param main an overall title for the plot
#' @param sub a sub title for the plot
#' @return Plot of mean value and corresponding probabilities for all simulated outcomes.
#' @details The default probabilty equals to 1/n. All the assigned probabilites must between 0 and 1.
#' @examples distr.simu.plot(distr = c(1,0.2,3.4,5,6.6,1.1,5,4.7,2.33,3), n = 4, times = 1000, col = 'red')
#' @export
#' @importFrom graphics curve hist par plot
#' @importFrom stats aggregate dnorm qqline qqnorm rbeta rbinom rchisq rexp rgamma rgeom rhyper rnbinom rnorm rpois runif var
#'

distr.simu.plot <- function (distr, n, times, prob = NULL, qqplot = FALSE, col = "black",
                        type = NULL, main = NULL, sub = NULL)
{
  if (is.null(prob)) {
    prob <- rep(1/length(distr), length(distr))
  }
  else {
    if (!identical(length(distr), length(prob))) {
      stop("'prob' is not the same length as 'possible outcomes'")
    }
    if (any(prob < 0)) {
      stop("'prob' contains negative values")
    }
    plist <- list()
    for (i in 1:n) {
      plist[[i]] = prob
    }
  }
  if (qqplot == FALSE){
    par(mfrow=c(1,1))
    RM <- rowMeans(replicate(n, sample(distr, size = times, replace = T,
                                       prob = prob)))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean")
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, mean=mean(distr), sd = sqrt(var(distr)/n)),col="red", add=TRUE)
  }
  if (qqplot == TRUE){
    par(mfrow=c(1,2))
    hist(RM,  main = "Histogram", freq = FALSE, xlab = "mean")
    x=seq(min(RM),max(RM),length=100)
    curve(dnorm(x, mean=mean(distr), sd = sqrt(var(distr)/n)),col="red", add=TRUE)
    qqnorm(RM)
    qqline(RM)
  }
}
