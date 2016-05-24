summary.lmf <-
function (object,
                         what.level = c("age", "year", "total"),
                         ...)
{
  #Rename the object
  z <- object
  #Collect call
  ans <- z[c("call")]
  #Choose "total" if all levels are choosen
  ifelse((all(what.level == c("age", "year", "total")) |
            all(what.level != c("age", "year", "total"))),
         ans$what.level <- "total",
         ans$what.level <- what.level)
  #Collect unique age classes and numbers of age classes
  ans$uage <- z$uage
  ans$nage <- z$nage
  #Collect unique years and numbers of years
  ans$uyear <- z$uyear
  ans$nyear <- z$nyear
  #Collect projection matrix
  ans$l <- z$l
  #Collect lambda, u and v
  ans$lambda <- z$lambda
  ans$u <- z$u
  ans$v <- z$v
  #Collect variance components with associated standard deviation and
  #degrees of freedom
  ans$sigma2.e <- z$sigma2.e
  ans$sigma2.dj <- z$sigma2.dj
  ans$sigma2.dj.sd <- z$sigma2.dj.sd
  ans$sigma2.dj.dof <- z$sigma2.dj.dof
  ans$sigma2.d <- z$sigma2.d
  ans$sigma2.d.sd <- z$sigma2.d.sd
  ans$sigma2.d.dof <- z$sigma2.d.dof
  #Check if estimates by age and year should be calculated
  if(ans$what.level == "age")
  {
    #Extract yearly alpha estimates by age with se, t-statistic and p-value
    est.ajt <- z$ajt
    se.ajt <- mapply(function(a) {lapply(a, function(b)sqrt(diag(b)))},
                     z$Ajt, SIMPLIFY = FALSE)
    tval.ajt <- mapply(function(a, b) {as.list(unlist(a)/unlist(b))}, est.ajt, se.ajt,
                       SIMPLIFY = FALSE)
    rdf.ajt <- mapply(as.list, z$dof, SIMPLIFY = FALSE)
    pval.ajt <- mapply(function(a, b){mapply(function(a, b)
    {2 * pt(abs(a), b, lower.tail = FALSE)}, a, b, SIMPLIFY = FALSE)},
                       tval.ajt, rdf.ajt, SIMPLIFY = FALSE)
    ans$coefficients.ajt <- cbind(unlist(est.ajt), unlist(se.ajt),
                                  unlist(tval.ajt), unlist(pval.ajt))
    rnames.ajt <- unlist(lapply(z$uage, function(a)
    {sapply(z$uyear, function(b) paste("age ", a, ", year ", b, ": ",
                                       names(z$ajt[[1]][[1]]), sep = ""))}))
    dimnames(ans$coefficients.ajt) <- list(rnames.ajt,
                                           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    #Extract yearly covariance matrices by age
    ans$Ajt <- z$Ajt
  }
  #Check if estimates by year should be calculated
  if(ans$what.level == "year" | ans$what.level == "age")
  {
    #Extract yearly covariance matrices
    ans$At <- lapply(z$At, function(a) a + z$M)
    #Extract yearly alpha estimates with se, t-statistic and p-value
    est.at <- z$at
    se.at <- lapply(ans$At, function(a)sqrt(diag(a)))
    tval.at <- as.list(unlist(est.at)/unlist(se.at))
    rdf.at <- as.list(Reduce('+', z$dof))
    pval.at <- mapply(function(a, b)
    {2 * pt(abs(a), b, lower.tail = FALSE)}, tval.at, rdf.at)
    ans$coefficients.at <- cbind(unlist(est.at), unlist(se.at),
                                 unlist(tval.at), pval.at)
    rnames.at <- sapply(z$uyear, function(a) paste("year ", a, ": ",
                                                   names(z$at[[1]]), sep = ""))
    dimnames(ans$coefficients.at) <- list(rnames.at,
                                          c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    #Extract yearly alpha estimates corrected for stampling error, with
    #se, t-statistic and p-value
    est.atC <- z$atC
    tval.atC <- as.list(unlist(est.atC)/unlist(se.at))
    pval.atC <- mapply(function(a, b)
    {2 * pt(abs(a), b, lower.tail = FALSE)}, tval.atC, rdf.at)
    ans$coefficients.atC <- cbind(unlist(est.atC), unlist(se.at),
                                  unlist(tval.atC), pval.atC)
    dimnames(ans$coefficients.atC) <- list(rnames.at,
                                           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  }
  #Check if total estimates should be calculated
  if(ans$what.level == "year" | ans$what.level == "age" | ans$what.level == "total")
  {
    #Extract temporal mean alpha estimates with se, t-statistic and p-value
    est.aM <- as.vector(z$aM)
    se.aM <-  as.vector(sqrt(diag(z$M)))
    tval.aM <- est.aM/se.aM
    rdf.aM <- z$sigma2.d.dof
    pval.aM <- 2 * pt(abs(tval.aM), rdf.aM, lower.tail = FALSE)
    ans$coefficients.aM <- cbind(est.aM, se.aM, tval.aM, pval.aM)
    dimnames(ans$coefficients.aM) <- list(dimnames(z$aM)[[2]],
                                          c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    #Extract temporal covariance matrix
    ans$M <- z$M
    #Extract alpha estimates assuming no fluctuating selection (M = 0),
    #with se, t-statistic and p-value
    est.anf <- as.vector(z$anf)
    se.anf <-  as.vector(sqrt(diag(z$Anf)))
    tval.anf <- est.anf/se.anf
    rdf.anf <- z$sigma2.d.dof
    pval.anf <- 2 * pt(abs(tval.anf), rdf.anf, lower.tail = FALSE)
    ans$coefficients.anf <- cbind(est.anf, se.anf, tval.anf, pval.anf)
    dimnames(ans$coefficients.anf) <- list(dimnames(z$anf)[[2]],
                                           c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    #Extract temporal covariance matrix
    ans$Anf <- z$Anf
  }
  #Define class
  class(ans) <- "summary.lmf"
  #Return
  ans
}
