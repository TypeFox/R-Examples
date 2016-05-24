summary.fitfrail <- function(object, type="survival", Lambda.times=NULL, 
                             censored=FALSE, se=FALSE, CI=0.95, ...) {
  fit <- object
  
  if (!match(type,c("survival", "hazard"), nomatch=0))
    stop("type must be either 'survival' or 'hazard'")
  
  if (!inherits(fit, "fitfrail")) 
    stop("summary.fitfrail can only be used for fitfrail objects")
  
  if (!is.null(Lambda.times)) {
    if (!is.numeric(Lambda.times)) 
      stop("Lambda.times must be numeric")
  } else if (censored) {
    Lambda.times <- fit$VARS$time
  } else {
    Lambda.times <- fit$VARS$time[fit$VARS$status > 0]
  }
  Lambda.times <- sort(unique(Lambda.times))
  
  stopifnot((CI > 0)&&(CI < 1))
  
  status <- as.integer(fit$VARS$status > 0)
  n.risk.total <- length(status)
  
  n.risk <- n.risk.total - vapply(Lambda.times, function(t) {
    sum(fit$VARS$time < t)
  }, 0) # num failures at time t-
  
  n.event <- c(sum(status[fit$VARS$time <= Lambda.times[1]]), 
               diff(vapply(Lambda.times, 
                 function(t) {
                   sum(status[fit$VARS$time <= t])
                 }, 0))) # num failures at time t-
  
  result <- data.frame(time=Lambda.times,
                       n.risk=n.risk,
                       n.event=n.event)
  rownames(result) <- NULL
  
  if (type == "survival") {
    result$surv <- exp(-fit$Lambda.fun(Lambda.times))
  } else if (type == "hazard") {
    result$haz <- fit$Lambda.fun(Lambda.times)
  }
  
  if (se) {
    cbh.se <- diag(vcov(fit, Lambda.times=Lambda.times, boot=TRUE, ...))
    cbh.se <- cbh.se[grepl("^Lambda", names(cbh.se))]
    
    if (type == "survival") {
      result$std.err <- exp(-result$surv + cbh.se^2/2)*sqrt(exp(cbh.se^2) - 1)
      
      if (CI > 0) {
        zval <- qnorm(1- (1-CI)/2, 0,1)
        lower <- pmax(result$surv - zval* result$std.err, 0)
        upper <- pmin(result$surv + zval* result$std.err, 1)
        result$lower.ci <- lower
        result$upper.ci <- upper
      }
    } else if (type == "hazard") {
      result$std.err <- cbh.se
      
      if (CI > 0) {
        zval <- qnorm(1- (1-CI)/2, 0,1)
        lower <- pmax(result$haz - zval* result$std.err, 0)
        upper <- result$haz + zval* result$std.err
        result$lower.ci <- lower
        result$upper.ci <- upper
      }
    }
    
  }
  
  result
}