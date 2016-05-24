AE3 <- function(model, random, data, seed=1234, n.sim=50000, verbose=TRUE)
{
  for(p in c("MASS", "nlme")) {
     if (length(grep(paste("^package:", p, "$", sep=""), search())) == 0) {
        if (!require(p, quietly = TRUE, character.only=TRUE))
        warning(paste("AE3 needs package `", p, "' to be fully functional; please install", sep=""))
     }
  }
  res <- nlme::lme(model, random = random, data = data, method = "ML")
  lns2 <- attr(res$apVar, "Pars")[1]
  lns1 <- attr(res$apVar, "Pars")[2]
  var22 <- res$apVar[1, 1]
  cov21 <- res$apVar[2, 1]
  var11 <- res$apVar[2, 2]
  h2 <- exp(2*lns2)/(exp(2*lns1) + exp(2*lns2))
  deriv1 <- -2*exp(2*lns1)*exp(2*lns2)/(exp(2*lns1) + exp(2*lns2))^2
  deriv2 <- -2*exp(4*lns2)/(exp(2*lns1) + exp(2*lns2))^2 + 2*exp(2*lns2)/(exp(2*lns1) + exp(2*lns2))
  se <- sqrt(var11*deriv1^2 + var22*deriv2^2 + 2*cov21*deriv1*deriv2)
  lcl95.no <- h2 - qnorm(.975)*se
  ucl95.no <- h2 + qnorm(.975)*se
  se.qnorm <- se/dnorm(qnorm(h2))
  lcl95.qnorm <- pnorm(qnorm(h2) - qnorm(.975)*se.qnorm)
  ucl95.qnorm <- pnorm(qnorm(h2) + qnorm(.975)*se.qnorm)
  diff <- lns1 - lns2
  se.diff <- sqrt(var11 + var22 - 2*cov21)
  diff.lower <- diff - 1.96*se.diff
  diff.upper <- diff + 1.96*se.diff
  h2.upper <- 1/(exp(2*diff.lower) + 1)
  h2.lower <- 1/(exp(2*diff.upper) + 1)
  set.seed(seed)
  mu = c(lns1, lns2)
  Sigma = matrix(c(var11, cov21, cov21, var22), ncol = 2)
  samp <- exp(2*MASS::mvrnorm(n.sim, mu, Sigma, empirical = TRUE))
  samp.h2 <- samp[,2]/(samp[,1] + samp[,2])
  CI <- quantile(samp.h2, c(0.025, 0.975))
  CL <- matrix(c(lcl95.no,ucl95.no, lcl95.qnorm,ucl95.qnorm, h2.lower,h2.upper, CI[1],CI[2]),ncol=2,byrow=TRUE)

  if(verbose)
  {
    print(summary(res))
    cat("\nHeritability:", h2,
        "\n95% CI (no transformation):     ", CL[1,1], "-", CL[1,2],
        "\n95% CI (probit transformation): ", CL[2,1], "-", CL[2,2],
        "\n95% CI (simpler transformation):", CL[3,1], "-", CL[3,2],
        "\n95% CI (simulated):             ", CL[4,1], "-", CL[4,2], "\n")
  }
  invisible(list(lme.result=res,h2=h2,CL=CL))
}
