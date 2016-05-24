  #################################################
  ## Fitting and Graphing Discrete Distributions ##
  #################################################

  data(HorseKicks)
  barplot(HorseKicks, col = 2,
          xlab = "Number of Deaths", ylab = "Number of Corps-Years",
          main = "Deaths by Horse Kicks")

  data(Federalist)
  barplot(Federalist, col = 2,
          xlab = "Occurrences of 'may'", ylab = "Number of Blocks of Text",
          main = "'may' in Federalist papers")

  data(WomenQueue)
  barplot(WomenQueue, col = 2,
          xlab = "Number of women", ylab = "Number of queues",
          main = "Women in queues of length 10")

  data(WeldonDice)
  barplot(WeldonDice, names = c(names(WeldonDice)[-11], "10+"), col = 2,
          xlab = "Number of 5s and 6s", ylab = "Frequency",
	  main = "Weldon's dice data")

  data(Butterfly)
  barplot(Butterfly, col = 2,
          xlab = "Number of individuals", ylab = "Number of Species",
  	main = "Butterfly species im Malaya")


  ############################
  ## Binomial distributions ##
  ############################

  par(mfrow = c(1,2))
  barplot(dbinom(0:10, p = 0.15, size = 10), names = 0:10, col = grey(0.7),
          main = "p = 0.15", ylim = c(0,0.35))
  barplot(dbinom(0:10, p = 0.35, size = 10), names = 0:10, col = grey(0.7),
          main = "p = 0.35", ylim = c(0,0.35))
  par(mfrow = c(1,1))
  mtext("Binomial distributions", line = 2, cex = 1.5)

  plot(0:10, dbinom(0:10, p = 0.15, size = 10), type = "b", ylab = "Density",
       ylim = c(0, 0.4), main = "Binomial distributions, N = 10", pch = 19)
  lines(0:10, dbinom(0:10, p = 0.35, size = 10), type = "b", col = 2, pch = 19)
  lines(0:10, dbinom(0:10, p = 0.55, size = 10), type = "b", col = 4, pch = 19)
  lines(0:10, dbinom(0:10, p = 0.75, size = 10), type = "b", col = 3, pch = 19)
  legend(3, 0.4, c("p", "0.15", "0.35", "0.55", "0.75"), lty = rep(1,5), col =
         c(0,1,2,4,3), bty = "n")

  ###########################
  ## Poisson distributions ##
  ###########################

  par(mfrow = c(1,2))
  dummy <- barplot(dpois(0:12, 2), names = 0:12, col = grey(0.7), ylim = c(0,0.3),
        main = expression(lambda == 2))
  abline(v = dummy[3], col = 2)
  diff <- (dummy[3] - dummy[2]) * sqrt(2)/2
  lines(c(dummy[3] - diff, dummy[3] + diff), c(0.3, 0.3), col = 2)
  dummy <- barplot(dpois(0:12, 5), names = 0:12, col = grey(0.7), ylim = c(0,0.3),
          main = expression(lambda == 5))
  abline(v = dummy[6], col = 2)
  diff <- (dummy[6] - dummy[5]) * sqrt(5)/2
  lines(c(dummy[6] - diff, dummy[6] + diff), c(0.3, 0.3), col = 2)
  par(mfrow = c(1,1))
  mtext("Poisson distributions", line = 2, cex = 1.5)

  #####################################
  ## Negative binomial distributions ##
  #####################################

  nbplot <- function(p = 0.2, size = 2, ylim = c(0, 0.2))
  {
    plot(0:20, dnbinom(0:20, p = p, size = size), type = "h", col = grey(0.7),
         xlab = "Number of failures (k)", ylab = "Density", ylim = ylim,
         yaxs = "i", bty = "L")
    nb.mean <- size * (1-p)/p
    nb.sd <- sqrt(nb.mean/p)
    abline(v = nb.mean, col = 2)
    lines(nb.mean + c(-nb.sd, nb.sd), c(0.01, 0.01), col = 2)
    legend(14, 0.2, c(paste("p = ", p), paste("n = ", size)), bty = "n")
  }
  par(mfrow = c(3,2))
  nbplot()
  nbplot(size = 4)
  nbplot(p = 0.3)
  nbplot(p = 0.3, size = 4)
  nbplot(p = 0.4, size = 2)
  nbplot(p = 0.4, size = 4)
  par(mfrow = c(1,1))
  mtext("Negative binomial distributions for the number of trials to observe n = 2 or n = 4 successes", line = 3)

  #####################
  ## Goodness of fit ##
  #####################

  p <- weighted.mean(as.numeric(names(HorseKicks)), HorseKicks)
  p.hat <- dpois(0:4, p)
  expected <- sum(HorseKicks) * p.hat
  chi2 <- sum((HorseKicks - expected)^2/expected)
  pchisq(chi2, df = 3, lower = FALSE)
  ## or:
  HK.fit <- goodfit(HorseKicks)
  summary(HK.fit)

  ## Are the dice fair?
  p.hyp <- 1/3
  p.hat <- dbinom(0:12, prob = p.hyp, size = 12)
  expected <- sum(WeldonDice) * p.hat
  expected <- c(expected[1:10], sum(expected[11:13]))
  chi2 <- sum((WeldonDice - expected)^2/expected)
  G2 <- 2*sum(WeldonDice*log(WeldonDice/expected))
  pchisq(chi2, df = 10, lower = FALSE)
  ## Are the data from a binomial distribution?
  p <- weighted.mean(as.numeric(names(WeldonDice))/12, WeldonDice)
  p.hat <- dbinom(0:12, prob = p, size = 12)
  expected <- sum(WeldonDice) * p.hat
  expected <- c(expected[1:10], sum(expected[11:13]))
  chi2 <- sum((WeldonDice - expected)^2/expected)
  G2 <- 2*sum(WeldonDice*log(WeldonDice/expected))
  pchisq(chi2, df = 9, lower = FALSE)
  ## or:
  WD.fit1 <- goodfit(WeldonDice, type = "binomial", par = list(prob = 1/3, size = 12))
  WD.fit1$fitted[11] <- sum(predict(WD.fit1, newcount = 10:12))
  WD.fit2 <- goodfit(WeldonDice, type = "binomial", par = list(size = 12), method = "MinChisq")
  summary(WD.fit1)
  summary(WD.fit2)

  F.fit1 <- goodfit(Federalist)
  F.fit2 <- goodfit(Federalist, type = "nbinomial")
  summary(F.fit1)
  par(mfrow = c(2,2))
  plot(F.fit1, scale = "raw", type = "standing")
  plot(F.fit1, type = "standing")
  plot(F.fit1)
  plot(F.fit1, type = "deviation")
  par(mfrow = c(1,1))
  plot(F.fit2, type = "deviation")
  summary(F.fit2)

  data(Saxony)
  S.fit <- goodfit(Saxony, type = "binomial", par = list(size = 12))
  summary(S.fit)
  plot(S.fit)

  ###############
  ## Ord plots ##
  ###############

  par(mfrow = c(2,2))
  Ord_plot(HorseKicks, main = "Death by horse kicks")
  Ord_plot(Federalist, main = "Instances of 'may' in Federalist papers")
  Ord_plot(Butterfly, main = "Butterfly species collected in Malaya")
  Ord_plot(WomenQueue, main = "Women in queues of length 10")
  par(mfrow = c(1,1))

  ###############
  ## Distplots ##
  ###############

  distplot(HorseKicks, type = "poisson")
  distplot(HorseKicks, type = "poisson", lambda = 0.61)
  distplot(Federalist, type = "poisson")
  distplot(Federalist, type = "nbinomial")
  distplot(Saxony, type = "binomial", size = 12)
