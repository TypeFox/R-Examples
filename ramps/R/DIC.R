DIC <- function(object, ...)
{
   UseMethod("DIC")
}

DIC.ramps <- function(object, ...)
{
  ## Posterior parameter means
   w <- rep(1 / nrow(object$params), nrow(object$params))
   params <- as.vector(crossprod(object$params, w))
   loglik <- as.numeric(crossprod(object$loglik, w))
   
   BETA <- params2beta(params, object$control)

   switch(object$control$mpdfun,
      mpdbeta = {
         val <- mpdbeta(params, object$y, object$xmat, object$kmat, object$wmat,
                        object$correlation, object$etype, object$ztype,
                        object$retype, object$weights, object$control)
      },
      mpdbetaz = {
         if(!all(object$control$z$monitor))
           stop("All z must be monitored to compute DIC")
         
         BETA <- c(BETA, as.vector(object$kmat %*% crossprod(object$z, w)))

         sites <- unique.sites(object$kmat)
         xk1mat <- cBind(object$xmat, sites$map)
         k2mat <- sites$coords

         val <- mpdbetaz(params, object$y, xk1mat, k2mat, object$wmat,
                         object$correlation, object$etype, object$ztype,
                         object$retype, object$weights, object$control)
      }
   )
   loglik.mu <-  -0.5 * length(object$y) * log(2.0 * pi) - val$logsqrtdet -
                    (crossprod(val$uXtSiginvX %*% (val$betahat - BETA))[1]
                     + val$quadform[1]) / 2.0

   ## Deviance information criterion
   pD <- -2.0 * (loglik - loglik.mu)
   c(DIC = pD - 2.0 * loglik, pD = pD)
}
