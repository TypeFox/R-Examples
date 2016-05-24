#' Pooled digital PCR analysis
#' 
#' Estimates mean number of template molecules per partition and concentration of sample 
#' from pooled replicates of experiments.
#' 
#' @aliases analyze_pooled
#' @param input object of class \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}}.
#' @param v volume (microliters).
#' @param conf.level confidence level of the intervals and groups.
#' @export
#' @note This function was implemented using the code in supplemental materials in 
#' Dorazio, 2015 (see References).
#' @return data frame with the number of rows equal to the number of replicates.
#' @author Robert M. Dorazio, Margaret E. Hunter.
#' @references Dorazio RM, Hunter ME, \emph{Statistical Models for the Analysis 
#' and Design of Digital Polymerase Chain Reaction (dPCR) Experiments}. 
#' Analytical Chemistry 2015. 87(21): p.10886-10893

analyze_pooled <- function(input, v, conf.level = 0.05) {
  # functionality and code below are taken from: 
  # Dorazio, R. M.; Hunter, M. E. Anal. Chem. 2015, 87 (21), 10886-10893.
  dat <- summary(input, print = FALSE)

  comp_data <- data.frame(exper = slot(input, "exper"), 
                          assay = slot(input, "assay"),
                          k = dat[["partitions"]][["k"]], 
                          n = dat[["partitions"]][["n"]])
  
  do.call(rbind, lapply(levels(comp_data[["exper"]]), function(single_experiment) {
    
    y <- comp_data[comp_data[["exper"]] == single_experiment, "k"]
    m <- comp_data[comp_data[["exper"]] == single_experiment, "n"]
    target <- as.character(comp_data[comp_data[["exper"]] == single_experiment, "assay"])
      
    v.offset <- rep(log(v), length(y))
    ymat <- cbind(y, m - y)
    fit <- glm(ymat ~ target - 1, family = binomial(link = "cloglog"), offset = v.offset)
    beta.mle <- fit[["coefficients"]]
    beta.vcv <- vcov(fit)
    zcrit <- qnorm(1 - conf.level/2)
    deviance <- fit[["deviance"]]

    Xvec <- matrix(c(-1, 1), ncol = 1)
    logR.est = t(Xvec) %*% beta.mle
    logR.var <- t(Xvec) %*% beta.vcv %*% Xvec
    
    data.frame(sample_name = as.character(single_experiment),
               R_est = exp(logR.est), 
               R_lowerCL = exp(logR.est - zcrit * sqrt(logR.var)),
               R_upperCL = exp(logR.est + zcrit * sqrt(logR.var)), 
               GOF = 1 - pchisq(deviance, df = length(y) - length(beta.mle)))
  }))
}




# negative of gradient of log-likelihood function

negGradLL <-  function(beta, y, m, v, X) {
  lambda <- as.vector(exp(X %*% beta))
  wvec <- lambda * (y/(1-exp(-lambda*v)) - m)
  retVal <- -v * as.vector(t(X) %*% wvec)
  retVal
}

# negative of log-likelihood function

negLL <- function(beta, y, m, v, X) {
  lambda <- as.vector(exp(X %*% beta))
  logL <- y*log(1-exp(-lambda*v)) - lambda*v*(m-y)
  (-1)*sum(logL)
}



# negative of hessian of log-likelihood function

negHessLL <-  function(beta, y, m, v, X) {
  lambda <- as.vector(exp(X %*% beta))
  uvec <- lambda * ( y*(1-(1+lambda*v)*exp(-lambda*v))/(1-exp(-lambda*v))^2 - m)
  retVal <- -v * t(X) %*% diag(uvec) %*% X 
  retVal
}

# negative of saturated log-likelihood function (used to compute deviance statistic)

negLL.saturated <- function(y, m, v) {
  lambda <- -(1/v)*log(1 - y/m)
  logL <- y*log(1-exp(-lambda*v)) - lambda*v*(m-y)
  (-1)*sum(logL)
}


