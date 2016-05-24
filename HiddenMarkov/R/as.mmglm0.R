"as.mmglm0" <-
function (object)
{
    #   turns dthmm object into mmpp
    x <- list(Pi=object$Pi, delta=object$delta,
              family=object$distn[2],
              link=object$distn[3],
              beta=rbind(object$pm$beta0, object$pm$beta1),
              sigma=object$pm$sigma,
              nonstat=object$nonstat)
    #   x$x <- list()
    #   x$x$y <- object$x
    #   x$x$x1 <- object$pn$x1
    #   x$x$size <- object$pn$size
    #   x$x <- as.data.frame(x$x)
    x$x <- object$glmdata
    x$glmformula <- object$glmformula
    class(x) <- c("mmglm0")
    return(x)
}

