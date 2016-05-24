ss.aipe.reliability <- function (model = NULL, type = NULL, width = NULL, S = NULL, 
            conf.level = 0.95, assurance = NULL, data = NULL, i = NULL, 
            cor.est = NULL, lambda = NULL, psi.square = NULL, initial.iter = 500, 
            final.iter = 5000, start.ss = NULL, verbose=FALSE) 
  {
    if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
    
    
    ############  
    #### This is an earlier ci.reliability() function (e.g., from version 3.2.0) that is used for this ss.aipe.reliability() function. 
    # The "ci.reliability" function originally used is now .original.ci.reliability() here so that ss.aipe.reliability() can use it. 
    # The new ci.reliability() function (used in Version 4.0.0) is more flexible and is changed too much to work with ss.aipe.reliability(). 
    
    .original.ci.reliability <- function (S = NULL, data = NULL, N = NULL, model = "True-Score Equivalent", 
                                          type = "Factor Analytic", conf.level = 0.95, interval = TRUE, Bootstrap = FALSE, B = 10000, BootstrapCI = "BCa") 
    {
      
      model.type1 <- c("Parallel", "SB", "Spearman Brown", "Spearman-Brown", 
                       "sb", "parallel")
      model.type2 <- c("True Score", "True Score Equivalent", "True-Score Equivalent", 
                       "Equivalent", "Tau Equivalent", "Chronbach", "Cronbach", 
                       "cronbach", "alpha", "true score", "tau-equivalent", 
                       "Tau-Equivalent", "True-Score", "true-score")
      model.type3 <- c("Congeneric", "congeneric", "omega", "Omega")
      if (sum(model == model.type1, model == model.type2, model == 
              model.type3) != 1) 
        stop("Assign one and only one of the three types of models to 'model'.")
      alpha <- 1 - conf.level
      if(Bootstrap == FALSE) {
        if (type == "Factor Analytic") {
          if (!is.null(data)) {
            S <- var(na.omit(data))
            N <- dim(data)[1]
          }
          if (!isSymmetric(S, tol = 1e-05)) 
            stop("Input a symmetric covariance or correlation matrix 'S.'")
          if (is.null(N)) 
            stop("Since only 'S' is entered, 'N' is also needed.")
          q <- nrow(S)
          if (sum(model == model.type1) == 1) 
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
          if (sum(model == model.type2) == 1) 
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
          if (sum(model == model.type3) == 1) 
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
          l <- length(result$Factor.Loadings)
          p <- length(result$Indicator.var)
          k <- nrow(result$Parameter.cov)
          if (l == 1) {
            u <- result$Factor.Loadings * q
            names(u) <- NULL
            var.u <- result$Parameter.cov[1, 1] * q * q
          }
          else {
            u <- sum(result$Factor.Loadings)
            var.u <- sum(result$Parameter.cov[1:q, 1:q])
          }
          if (p == 1) {
            v <- result$Indicator.var * q
            names(v) <- NULL
            var.v <- result$Parameter.cov[k, k] * q * q
          }
          else {
            v <- sum(result$Indicator.var)
            var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
          }
          if (l != 1 & p != 1) 
            cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
          if (l == 1 & p != 1) 
            cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
          if (l == 1 & p == 1) 
            cov.u.v <- result$Parameter.cov[2, 1] * q * q
          omega <- u^2/(u^2 + v)
          if(interval == "T" | interval == "t") {
            D1 <- 2 * u * v/(u^2 + v)^2
            D2 <- (-1) * u^2/(u^2 + v)^2
            se.omega <- sqrt(D1^2 * var.u + D2^2 * var.v + 2 * D1 * D2 * cov.u.v)
            z <- qnorm(1 - alpha/2)
            ci.lower <- omega - z * se.omega
            if (ci.lower < 0) {
              ci.lower = 0
            }
            ci.upper <- omega + z * se.omega
            if (ci.upper > 1) {
              ci.upper = 1
            }
          } else {
            ci.lower <- ci.upper <- se.omega <- -1
          }
          ci <- list(CI.lower = ci.lower, CI.upper = ci.upper, 
                     Estimated.reliability = omega, SE.reliability = se.omega, 
                     Conf.Level = conf.level)
          return(ci)
        }
        if (type == "Normal Theory") {
          if (sum(model == model.type3) == 1) 
            stop("The Congeneric model can not be used with `Normal Theory.'")
          if (!is.null(data)) { 
            data <- na.omit(data)
            if (dim(data)[1] == dim(data)[2]) {
              if (sum(round(data, 4) == t(round(data, 4))) == (dim(data)[1] * dim(data)[2])) 
                cov.data <- data
            }
            if (dim(data)[1] != dim(data)[2]) {
              cov.data <- var(data, y = NULL, na.rm = TRUE)
              N <- dim(data)[1]
            }
          }
          if (!is.null(S)) 
            cov.data <- S
          if (is.null(N)) 
            stop("Since only 'S' is entered, 'N' is also needed.")
          if (sum(model == model.type1) == 1) {
            sigma.jj <- sum(diag(cov.data))
            sigma2.Y <- sum(cov.data)
            p <- ncol(cov.data)
            alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
            if(interval) {
              Crit.Value <- qnorm(1 - (1 - conf.level)/2)
              variance <- (2 * (1 - alpha)^2 * p)/((N - 1) * (p - 1))
              SE <- sqrt((2 * (1 - alpha)^2 * p)/((N - 1) * (p - 1)))
              UB <- alpha + (Crit.Value * SE)
              if (UB > 1) {
                UB = 1
              }
              LB <- alpha - (Crit.Value * SE)
              if (LB < 0) {
                LB = 0
              }
            } else {
              LB <- UB <- SE <- -1
            }
            CI <- list(CI.lower = LB, CI.upper = UB, Estimated.reliability = alpha, 
                       SE.reliability = SE, Conf.Level = conf.level)
          }
          if (sum(model == model.type2) == 1) {
            sigma.jj <- sum(diag(cov.data))
            sigma2.Y <- sum(cov.data)
            p <- ncol(cov.data)
            alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
            if(interval == "T" | interval == "t") {
              cor.mat <- cov2cor(cov.data)
              p <- ncol(cor.mat)
              j <- cbind(rep(1, times = p))
              step.1 <- (p^2/(p - 1)^2)
              gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
              gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
              gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
              gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
              gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% cor.mat) %*% j)
              gamma.2 <- gamma.2.1 - gamma.2.2
              gamma.final <- gamma.1 * gamma.2
              n <- N - 1
              variance <- (step.1 * gamma.final)/n
              SE <- sqrt(variance)
              Crit.Value <- qnorm(1 - (1 - conf.level)/2)
              UB <- alpha + (Crit.Value * SE)
              if (UB > 1) {
                UB = 1
              }
              LB <- alpha - (Crit.Value * SE)
              if (LB < 0) {
                LB = 0
              }
            } else {
              LB <- UB <- SE <- -1
            }
            CI <- list(CI.lower = LB, CI.upper = UB, Estimated.reliability = alpha, 
                       SE.reliability = SE, Conf.Level = conf.level)
          }
          return(CI)
        }
      } else {
        if(!requireNamespace("boot", quietly = TRUE)) stop("The package 'boot' is needed; please install the package and try again.")
        BootstrapCI.type1 <- c("perc", "Perc", "percentile", "Percentile", "percentile ci", "Percentile ci", "Percentile CI", "percentile CI")
        BootstrapCI.type2 <- c("bca", "Bca", "BCa", "BCA", "Bias-corrected and acceleration", "Bias-Corrected and Acceleration", 
                               "bias-corrected and acceleration") 
        if (sum(BootstrapCI == BootstrapCI.type1, BootstrapCI == BootstrapCI.type2) != 1) 
          stop("Assign one and only one of the two types of percentile methodmodels to 'model'.")
        if (is.null(data))
          stop("Data is required for running bootstrap CI")
        if (type == "Factor Analytic") {
          model <- 0
          S <- var(data, y = NULL, na.rm = TRUE)
          N <- dim(data)[1]
          q <- nrow(S)
          if (sum(model == model.type1) == 1) {
            model <- 1
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
          } else if(sum(model == model.type2) == 1)  {
            model <- 2
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
          } else {
            model <- 3
            result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
          }
          l <- length(result$Factor.Loadings)
          p <- length(result$Indicator.var)
          k <- nrow(result$Parameter.cov)
          if (l == 1) {
            u <- result$Factor.Loadings * q
            names(u) <- NULL
            var.u <- result$Parameter.cov[1, 1] * q * q
          } else {
            u <- sum(result$Factor.Loadings)
            var.u <- sum(result$Parameter.cov[1:q, 1:q])
          }
          if (p == 1) {
            v <- result$Indicator.var * q
            names(v) <- NULL
            var.v <- result$Parameter.cov[k, k] * q * q
          } else {
            v <- sum(result$Indicator.var)
            var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
          }
          if (l != 1 & p != 1) 
            cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
          if (l == 1 & p != 1) 
            cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
          if (l == 1 & p == 1) 
            cov.u.v <- result$Parameter.cov[2, 1] * q * q
          omega <- u^2/(u^2 + v)
          .bs.omega <- function(data, model, i) {
            S <- var(data[i, ], y = NULL, na.rm = TRUE)
            N <- dim(data)[1]
            q <- nrow(S)
            if (model == 1) {
              result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = TRUE)
            } else if (model == 2) {
              result <- CFA.1(S = S, N = as.numeric(N), equal.loading = TRUE, equal.error = FALSE)
            } else {
              result <- CFA.1(S = S, N = as.numeric(N), equal.loading = FALSE, equal.error = FALSE)
            }
            l <- length(result$Factor.Loadings)
            p <- length(result$Indicator.var)
            k <- nrow(result$Parameter.cov)
            if (l == 1) {
              u <- result$Factor.Loadings * q
              names(u) <- NULL
              var.u <- result$Parameter.cov[1, 1] * q * q
            } else {
              u <- sum(result$Factor.Loadings)
              var.u <- sum(result$Parameter.cov[1:q, 1:q])
            }
            if (p == 1) {
              v <- result$Indicator.var * q
              names(v) <- NULL
              var.v <- result$Parameter.cov[k, k] * q * q
            } else {
              v <- sum(result$Indicator.var)
              var.v <- sum(result$Parameter.cov[(k - q + 1):k, (k - q + 1):k])
            }
            if (l != 1 & p != 1) 
              cov.u.v <- (sum(result$Parameter.cov) - var.v - var.u)/2
            if (l == 1 & p != 1) 
              cov.u.v <- sum(result$Parameter.cov[2:k, 1]) * q
            if (l == 1 & p == 1) 
              cov.u.v <- result$Parameter.cov[2, 1] * q * q
            omega <- u^2/(u^2 + v)    
            return(omega)
          }
          boot.out <- boot::boot(data = data, statistic = .bs.omega, R = B, stype = "i", model = model)
          se.omega <- sd(boot.out$t)
          CI.output <- NULL
          CI.lower <- NULL
          CI.upper <- NULL
          if (sum(BootstrapCI == BootstrapCI.type1) == 1) {
            CI.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
            CI.lower <- CI.output[4]
            CI.upper <- CI.output[5]
          } else { 
            CI.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
            CI.lower <- CI.output[4]
            CI.upper <- CI.output[5]
          }
          ci <- list(CI.lower = CI.lower, CI.upper = CI.upper, 
                     Estimated.reliability = omega, SE.reliability = se.omega, 
                     Conf.Level = conf.level)	 	
          return(ci)
        }
        if (type == "Normal Theory") {
          if (sum(model == model.type3) == 1) 
            stop("The Congeneric model can not be used with `Normal Theory.'")
          cov.data <- var(data, y = NULL, na.rm = TRUE)
          N <- dim(data)[1]
          sigma.jj <- sum(diag(cov.data))
          sigma2.Y <- sum(cov.data)
          p <- ncol(cov.data)
          alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
          .bs.alpha <- function(data, i) {
            cov.data <- var(data[i, ], y = NULL, na.rm = TRUE)
            N <- dim(data)[1]
            sigma.jj <- sum(diag(cov.data))
            sigma2.Y <- sum(cov.data)
            p <- ncol(cov.data)
            alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
            return(alpha)
          }
          boot.out <- boot::boot(data = data, statistic = .bs.alpha, R = B, stype = "i")
          se.alpha <- sd(boot.out$t)
          CI.output <- NULL
          CI.lower <- NULL
          CI.upper <- NULL
          if (sum(BootstrapCI == BootstrapCI.type1) == 1) {
            CI.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "perc")$perc
            CI.lower <- CI.output[4]
            CI.upper <- CI.output[5]
          } else { 
            CI.output <- boot::boot.ci(boot.out = boot.out, conf = conf.level, type = "bca")$bca
            CI.lower <- CI.output[4]
            CI.upper <- CI.output[5]
          }
          ci <- list(CI.lower = CI.lower, CI.upper = CI.upper, 
                     Estimated.reliability = alpha, SE.reliability = se.alpha, 
                     Conf.Level = conf.level)	 	
          return(ci)
        }
      }
    }
    ############  
    
    
    
    model.type1 <- c("Parallel", "SB", "Spearman Brown", "Spearman-Brown", 
                     "sb", "parallel")
    model.type2 <- c("True Score", "True Score Equivalent", "True-Score Equivalent", 
                     "Equivalent", "Tau Equivalent", "Chronbach", "Cronbach", 
                     "cronbach", "alpha", "true score", "tau-equivalent", 
                     "Tau-Equivalent", "True-Score", "true-score")
    model.type3 <- c("Congeneric", "congeneric", "omega", "Omega")
    if (sum(model == model.type1, model == model.type2, model == 
            model.type3) != 1) 
      stop("Assign one and only one of the three types of models to 'model'.")
    if (sum(model == model.type1) == 1) 
      Model.to.Use <- "Parallel"
    if (sum(model == model.type2) == 1) 
      Model.to.Use <- "True Score"
    if (sum(model == model.type3) == 1) 
      Model.to.Use <- "Congeneric"
    type.1 <- c("Normal Theory", "Normal theory", "normal theory", 
                "nt", "NT")
    type.2 <- c("Factor Analytic", "factor analytic", "Factor analytic", 
                "factor Analytic")
    if (sum(type == type.1, type == type.2) != 1) 
      stop("Assign either Factor Analytic or Normal Theory to 'type'.")
    if (sum(type == type.1) == 1) 
      Type.to.Use <- "Normal Theory"
    if (sum(type == type.2) == 1) 
      Type.to.Use <- "Factor Analytic"
    if (sum(model == model.type1) == 1) {
      if (!is.null(cor.est)) {
        if (!is.null(lambda)) 
          stop("Please enter either cor.est or lambda, but not both.")
        if (is.null(psi.square)) 
          stop("Problem: please enter all of the necessary information: `i', `cor.est' or `lambda', and `psi.square.'")
        if (psi.square <= 0) 
          stop("Problem: `psi.square' must be greater than zero")
        if (length(psi.square) > 1) 
          stop("Problem: 'psi.square' must be one number for the Parallel model. If you want to enter multiple 'psi.square' values, you must use either `True Score' or `Congeneric.'")
        if (is.null(i)) 
          stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
        if (i <= 0) 
          stop("Problem: `i' must be greater than zero")
        if (length(cor.est) >= 2) 
          stop("Problem: you can only enter one 'cor.est' value.")
        lambda.1 <- sqrt(cor.est)
        v.lambda <- matrix(data = lambda.1, nrow = 1, ncol = i)
        v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
        temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                                    tol.det = 1e-05)
        cov.mat <- Sigma <- temp.mat$Population.Covariance
      }
      if (!is.null(lambda)) {
        if (!is.null(cor.est)) 
          stop("Problem: please enter either cor.est or lambda, but not both.")
        if (is.null(psi.square)) 
          stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
        if (psi.square <= 0) 
          stop("Problem: `psi.square' must be greater than zero.")
        if (length(psi.square) > 1) 
          stop("Problem: 'psi.square' must be one number for the Parallel model. If you want to enter multiple 'psi.square' values, you must use either the True Score or Congeneric model.")
        if (is.null(i)) 
          stop("Problem: please enter all of the necessary information: i, cor.est or lambda, and psi.square.")
        if (i <= 0) 
          stop("Problem: `i' must be greater than zero.")
        v.lambda <- matrix(data = lambda, nrow = 1, ncol = i)
        v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
        temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                                    tol.det = 1e-05)
        cov.mat <- Sigma <- temp.mat$Population.Covariance
      }
      if (!is.null(S)) {
        if (!isSymmetric(S, tol = 1e-05)) 
          stop("Input a symmetric covariance or correlation matrix, 'S'")
        cov.mat <- S
        i <- dim(S)[1]
      }
      if (!is.null(data)) {
        data <- na.omit(data)
        cov.mat <- var(data, y = NULL, na.rm = TRUE)
      }
      sigma.jj <- sum(diag(cov.mat))
      sigma2.Y <- sum(cov.mat)
      p <- ncol(cov.mat)
      alpha <- (p/(p - 1)) * (1 - sigma.jj/sigma2.Y)
      k <- 1 - alpha
      Crit.Value <- qnorm(1 - (1 - conf.level)/2)
      Nec.N <- as.numeric(((Crit.Value^2) * 8 * (k^2) * p)/((width^2) * 
                                                              (p - 1)) + 1)
    }
    if (sum(model == model.type2) == 1) {
      if (!is.null(cor.est)) {
        if (!is.null(lambda)) 
          stop("Problem: please enter either `cor.est' or `lambda', but not both")
        if (is.null(psi.square)) 
          stop("Problem: please enter a vector of `psi.square' values")
        if (is.null(i)) 
          stop("Problem: please enter a value for i")
        if ((i) != length(psi.square)) 
          stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
        if (i <= 0) 
          stop("Problem: `i' must be greater than zero")
        if (length(cor.est) >= 2) 
          stop("You can only enter one 'cor.est' value.")
        lambda.1 <- sqrt(cor.est)
        lambda.vector <- rep(lambda.1, times = i)
        Population.Cov <- covmat.from.cfm(Lambda = lambda.vector, 
                                          Psi.Square = psi.square)$Population.Covariance
        cor.mat <- cov2cor(Population.Cov)
        cov.mat <- cor.mat
      }
      if (!is.null(lambda)) {
        if (!is.null(cor.est)) 
          stop("Please enter either cor.est or lambda, but not both")
        if (is.null(psi.square)) 
          stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
        if (length(psi.square) != (i)) 
          stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
        if (is.null(i)) 
          stop("You need to enter a quantity for i if you enter lambda and psi.square")
        if (i <= 0) 
          stop("i must be greater than zero")
        if (length(lambda) > 1) 
          stop("'lambda' must be one number for the True Score model. If you want to enter multiple 'lambda' values, you must use the Congeneric model.")
        lambda.vector <- rep(lambda, times = i)
        Population.Cov <- covmat.from.cfm(Lambda = lambda.vector, 
                                          Psi.Square = psi.square)$Population.Covariance
        cor.mat <- cov2cor(Population.Cov)
        cov.mat <- cor.mat
      }
      if (!is.null(S)) {
        if (!isSymmetric(S, tol = 1e-05)) 
          stop("Input a symmetric covariance or correlation matrix 'S'")
        cor.mat <- cov2cor(S)
        cov.mat <- cor.mat
        i <- dim(S)[1]
      }
      if (!is.null(data)) {
        data <- na.omit(data)
        cor.mat <- cor(data)
        cov.mat <- cor.mat
      }
      p <- ncol(cor.mat)
      j <- cbind(rep(1, times = p))
      Crit.Value <- qnorm(1 - (1 - conf.level)/2)
      step.1 <- (p^2/(p - 1)^2)
      gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
      gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
      gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
      gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
      gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% 
                                                           cor.mat) %*% j)
      gamma.2 <- gamma.2.1 - gamma.2.2
      gamma.final <- gamma.1 * gamma.2
      Nec.N <- as.numeric(((step.1 * gamma.final)/((width/(2 * 
                                                             Crit.Value))^2)) + 1)
    }
    if (sum(model == model.type3) == 1) {
      if (!is.null(cor.est)) {
        if (!is.null(lambda)) 
          stop("Please enter either cor.est or lambda, but not both")
        if (is.null(psi.square)) 
          stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
        if (is.null(i)) 
          stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
        if (i <= 0) 
          stop("i must be greater than zero")
        if ((i) != length(psi.square)) 
          stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
        if (length(cor.est) >= 2) 
          stop("If you have multiple values for 'cor.est', please put them as 'lambda' values. The square root of cor.est equals lambda.")
        print("You entered one value for 'cor.est' with the Congeneric model. This model allows for multiple 'lambda' values. If you have only one 'cor.est' or 'lambda' value, you might want to use the True Score model.")
        lambda.1 <- sqrt(cor.est)
        v.lambda <- matrix(data = lambda.1, nrow = 1, ncol = i)
        v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
        temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                                    tol.det = 1e-05)
        cov.mat <- Sigma <- temp.mat$Population.Covariance
        cor.mat <- cov2cor(cov.mat)
      }
      if (!is.null(lambda)) {
        if (!is.null(cor.est)) 
          stop("Please enter either cor.est or lambda, but not both")
        if (is.null(psi.square)) 
          stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
        if (is.null(i)) 
          stop("Please enter all of the necessary information: i, cor.est or lambda, and psi.square")
        if (i <= 0) 
          stop("i must be greater than zero")
        if ((i) != length(lambda)) 
          stop("The number of values entered in the lambda vector should be the same as the quantity entered for i")
        if ((i) != length(psi.square)) 
          stop("The number of values entered in the psi.square vector should be the same as the quantity entered for i")
        v.lambda <- matrix(data = lambda, nrow = 1, ncol = i)
        v.Psi <- matrix(data = psi.square, nrow = 1, ncol = i)
        temp.mat <- covmat.from.cfm(Lambda = v.lambda, Psi.Square = v.Psi, 
                                    tol.det = 1e-05)
        cov.mat <- Sigma <- temp.mat$Population.Covariance
        cor.mat <- cov2cor(cov.mat)
      }
      if (!is.null(data)) {
        data <- na.omit(data)
        cov.mat <- var(data, y = NULL, na.rm = TRUE)
        cor.mat <- cov2cor(cov.mat)
      }
      if (!is.null(S)) {
        if (!isSymmetric(S, tol = 1e-05)) 
          stop("Input a symmetric covariance or correlation matrix 'S'")
        cor.mat <- cov2cor(S)
        cov.mat <- cor.mat
        i <- dim(S)[1]
      }
      p <- ncol(cor.mat)
      j <- cbind(rep(1, times = p))
      Crit.Value <- qnorm(1 - (1 - conf.level)/2)
      step.1 <- (p^2/(p - 1)^2)
      gamma.1 <- 2/((t(j) %*% cor.mat %*% j)^3)
      gamma.2.1.1 <- (t(j) %*% cor.mat %*% j)
      gamma.2.1.2 <- ((sum(diag(cor.mat %*% cor.mat))) + (sum(diag(cor.mat)))^2)
      gamma.2.1 <- gamma.2.1.1 * gamma.2.1.2
      gamma.2.2 <- 2 * (sum(diag(cor.mat))) * (t(j) %*% (cor.mat %*% 
                                                           cor.mat) %*% j)
      gamma.2 <- gamma.2.1 - gamma.2.2
      gamma.final <- gamma.1 * gamma.2
      Nec.N <- as.numeric(((step.1 * gamma.final)/((width/(2 * 
                                                             Crit.Value))^2)) + 1)
    }
    if (!is.null(assurance)) {
      if (assurance > 1) 
        assurance <- assurance/100
      print("An a priori Monte Carlo simulation study has been started so that the exact sample size for the requested condition can be determined. Please be patient, as this process may take several minutes (or longer given the computer and condition).")
      initial.assurance.N <- ceiling(Nec.N) + 1
      if (sum(model == model.type1) == 1) 
        Model.to.Use <- "Parallel"
      if (sum(model == model.type2) == 1) 
        Model.to.Use <- "True Score"
      if (sum(model == model.type3) == 1) 
        Model.to.Use <- "Congeneric"
      n.i <- initial.assurance.N
      if (is.null(start.ss)) {
        Difference <- -1
        while (Difference < 0) {
          CI.Result <- rep(0, initial.iter)
          for (iters in 1:initial.iter) {
            sim.data <- var(MASS::mvrnorm(n = n.i, mu = rep(0, 
                                                            i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
            if (sum(model == model.type1) == 1) 
              CI.Result.raw <- try(.original.ci.reliability(S = sim.data, 
                                                            N = n.i, model = "Parallel", type = Type.to.Use, 
                                                            conf.level = conf.level), silent = TRUE)
            if (sum(model == model.type2) == 1) 
              CI.Result.raw <- try(.original.ci.reliability(S = sim.data, 
                                                            N = n.i, model = "True Score", type = Type.to.Use, 
                                                            conf.level = conf.level), silent = TRUE)
            if (sum(model == model.type3) == 1) 
              CI.Result.raw <- try(.original.ci.reliability(S = sim.data, 
                                                            N = n.i, model = "Congeneric", type = Type.to.Use, 
                                                            conf.level = conf.level), silent = TRUE)
            if (class(CI.Result.raw) == "try-error") {
              CI.Result[iters] <- NA
            }
            if (class(CI.Result.raw) != "try-error") {
              CI.Result[iters] <- CI.Result.raw$CI.upper - 
                CI.Result.raw$CI.lower
            }
            Difference <- mean(na.omit(CI.Result) < width) - assurance
            if(verbose==TRUE) print(paste("The current assurance is", mean(na.omit(CI.Result) < width), "at the current sample size of", n.i, "you are in the 'initial.iter' stage."))
            if (Difference < 0) break
          }
          if (Difference < 0) {
            n.i <- n.i + 1
          }
        }
        initial.n.i <- n.i
      }
      if (!is.null(start.ss)) {
        CI.Result <- rep(0,final.iter)
        n.i <- start.ss
        Difference <- 0
        for (iters in 1:final.iter) {
          sim.data <<- var(MASS::mvrnorm(n = n.i, mu = rep(0, 
                                                           i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
          CI.Result.raw <<- try(.original.ci.reliability(S = sim.data, 
                                                         N = n.i, model = Model.to.Use, type = Type.to.Use, 
                                                         conf.level = conf.level), silent = TRUE)
          if (class(CI.Result.raw) == "try-error") {
            CI.Result[iters] <- NA
          }
          if (class(CI.Result.raw) != "try-error") {
            CI.Result[iters] <- CI.Result.raw$CI.upper - 
              CI.Result.raw$CI.lower
          }
          Difference <- mean(na.omit(CI.Result) < width) - assurance
          if(verbose==TRUE) print(paste("The current assurance is", mean(na.omit(CI.Result) < width), "at the current sample size of", n.i, "you are in the 'final.iter' stage."))
          if (Difference < 0) break
        }
        
        if (Difference > 0) 
          while (Difference > 0) {
            CI.Result <- rep(1, final.iter)
            for (iters in 1:final.iter) {
              sim.data <- var(MASS::mvrnorm(n = n.i, mu = rep(0, 
                                                              i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
              CI.Result.raw <- try(.original.ci.reliability(S = sim.data, 
                                                            N = n.i, model = Model.to.Use, type = Type.to.Use, 
                                                            conf.level = conf.level), silent = TRUE)
              if (class(CI.Result.raw) == "try-error") {
                CI.Result[iters] <- NA
              }
              if (class(CI.Result.raw) != "try-error") {
                CI.Result[iters] <- CI.Result.raw$CI.upper - 
                  CI.Result.raw$CI.lower
              }
              Difference <- mean(na.omit(CI.Result) < width) - assurance
              if(verbose==TRUE) print(paste("The current assurance is", mean(na.omit(CI.Result) < width), "at the current sample size of", n.i, "you are in the 'final.iter' stage."))
              if (Difference > 0) break
            }
            
            if (Difference > 0) {
              n.i <- n.i - 1
            }
          }
      }
      Difference <- -1
      while (Difference < 0) {
        CI.Result <- rep(0, final.iter)
        for (iters in 1:final.iter) {
          sim.data <- var(MASS::mvrnorm(n = n.i, mu = rep(0, 
                                                          i), Sigma = cov.mat, tol = 1e-06, empirical = FALSE))
          CI.Result.raw <- try(.original.ci.reliability(S = sim.data, 
                                                        N = n.i, model = Model.to.Use, type = Type.to.Use, 
                                                        conf.level = conf.level), silent = TRUE)
          if (class(CI.Result.raw) == "try-error") {
            CI.Result[iters] <- NA
          }
          if (class(CI.Result.raw) != "try-error") 
            {
            CI.Result[iters] <- CI.Result.raw$CI.upper - 
              CI.Result.raw$CI.lower
          }
          Difference <- mean(na.omit(CI.Result) < width) - assurance
          if(verbose==TRUE) print(paste("The current assurance is", mean(na.omit(CI.Result) < width), "at the current sample size of", n.i, "; you are in the 'final.iter' stage."))
          if (Difference < 0) break
        }
        
        if (Difference < 0) {
          n.i <- n.i + 1
        }
      }
      Nec.N.assurance <- n.i
      empirical.assurance <- round(mean(na.omit(CI.Result) < 
                                          width), 3)
      print(paste("A sample size of", n.i, "leads to an empirical assurance of", 
                  round(mean(na.omit(CI.Result) < width), 3)))
    }
    if (is.null(assurance)) 
      return(list(Required.Sample.Size = ceiling(Nec.N)))
    if (!is.null(assurance)) 
      return(list(Required.Sample.Size = ceiling(Nec.N.assurance), 
                  width = width, specified.assurance = assurance, empirical.assurance = empirical.assurance, 
                  final.iter = final.iter))
  }