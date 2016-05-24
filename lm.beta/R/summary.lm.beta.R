summary.lm.beta <- function(object, standardized=TRUE, ...) {
  object2 <- object
  object2$standardized.coefficients <- NULL
  attr(object2,"class") <- "lm"
  x.summary <- summary(object2,...)
  if(standardized) {
    x.summary$coefficients <- cbind(x.summary$coefficients[,1,drop=F],Standardized=object$standardized.coefficients[rownames(x.summary$coefficients)],x.summary$coefficients[,-1,drop=F])
    attr(x.summary,"class") <- c("summary.lm.beta","summary.lm")
  }
  x.summary
}