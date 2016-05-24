summary.CADFtest <- function(object, ...)
{
  # object is an object of class CADFtest
  rnames <- 
    c("t-test statistic:                         ",
      "estimated rho^2:                          ",
      "p-value:                                  ",
      "Max lag of the diff. dependent variable:  ",
      "Max lag  of the stationary covariate(s):  ",
      "Max lead of the stationary covariate(s):  ")

  cnames <- "CADF test"

  if (is.null(object$parameter))
  {
    rnames <- rnames[c(1,3:4)]
    cnames <- "ADF test"
  }

  test.summary <- matrix(NA,(6-3*as.numeric(is.null(object$parameter))), 1, 
	dimnames=list(rnames,cnames))

  test.summary[1] <- object$statistic
  test.summary[3-as.numeric(is.null(object$parameter))] <- object$p.value
  test.summary[4-as.numeric(is.null(object$parameter))] <- object$max.lag.y

  if (!is.null(object$parameter))
  {
    test.summary[2] <- object$parameter
    test.summary[5] <- object$max.lag.X
    test.summary[6] <- object$min.lag.X
  }
  
    model.summary <- summary.lm(object$est.model)
    pos.lag.dep <- which(rownames(model.summary$coefficients)=="L(y, 1)")
    model.summary$coefficients[pos.lag.dep, 4] <- object$p.value

    F <- NA
    df.num <- NA
    df.den <- NA
    k1 <- 0; k0 <- 0

    k <- dim(object$est.model$model)[2]-1
    T <- dim(object$est.model$model)[1]

    if ((object$type == "trend") & (k > 3))
    {
      reduced.model <- lm(object$est.model$model[,1] ~ object$est.model$model[,2] + object$est.model$model[,3])
      k1 <- k + 1; k0 <- 3
    }

    if ((object$type=="drift") & (k > 2))
    {
      reduced.model <- lm(object$est.model$model[,1] ~ object$est.model$model[,2])
      k1 <- k + 1; k0 <- 2
    }

    if ((object$type=="none") & (k > 2))
    {
      reduced.model <- lm(object$est.model$model[,1] ~ -1 + object$est.model$model[,2])
      k1 <- k; k0 <- 1 
    }

    if (k1 > k0)
    {
      s.reduced <- summary(reduced.model)

      RSS1 <- sum(model.summary$residuals^2)
      RSS0 <- sum(s.reduced$residuals^2)

      df.num <- k1 - k0
      df.den <- T - k1

      F = ((RSS0 - RSS1)/df.num) / (RSS1/df.den)
    }

    model.summary$fstatistic[1:3] <- c(F, df.num, df.den)

    CADFtestsummary <- list(test.summary=test.summary,
                            model.summary=model.summary)

    class(CADFtestsummary) <- c("CADFtestsummary", "summary.dynlm", "summary.lm")  
    return(CADFtestsummary)
}
