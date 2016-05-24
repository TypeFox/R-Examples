interval <- function(glm.object, ...)
  UseMethod("interval")

interval.glm <- function(glm.object, linkfit.object,
                         type = c("link", "response"),
                         conf.level = 0.95, ...) {
  if (missing(linkfit.object))
    linkfit.object <- predict(glm.object, type="link", se.fit=TRUE, ...)
  type <- match.arg(type)
  crit.value <- qt(1-(1-conf.level)/2, glm.object$df.residual)

  conf.int <- {
    half.conf.interval <- crit.value * linkfit.object$se.fit
    cbind(ci.low=linkfit.object$fit - half.conf.interval,
          ci.hi =linkfit.object$fit + half.conf.interval)
  }

  pred.int <- {
    half.pred.interval <- crit.value *
      sqrt(linkfit.object$se.fit^2 +
           linkfit.object$residual.scale^2)
    cbind(pi.low=linkfit.object$fit - half.pred.interval,
          pi.hi =linkfit.object$fit + half.pred.interval)
  }
  result <- cbind(fit=linkfit.object$fit, conf.int, pred.int)
  if (type=="link") result
  else family(glm.object)$linkinv(result)
}
