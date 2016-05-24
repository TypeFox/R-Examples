Rsq <-
function(model)
{
   if (!inherits(model, c("lm", "aov", "nls")))
      stop ("'Rsq' is only applied to the classes: 'lm', 'aov' or 'nls'.")
   if (inherits(model, c("glm", "manova", "maov", "mlm")))
      stop("'Rsq' is not applied to an object of this class!")

   pred <- predict(model)
   n <- length(pred)
   res <- resid(model)
   w <- weights(model)
   if (is.null(w)) w <- rep(1, n)
   rss <- sum(w * res ^ 2)
   resp <- pred + res
   center <- weighted.mean(resp, w)
   if (inherits(model, c("lm", "aov"))) {
      r.df <- model$df.residual
      int.df <- attr(model$terms, "intercept")
      if (int.df) {
         mss <- sum(w * scale(pred, scale = FALSE)^2)
      } else {
         mss <- sum(w * scale(pred, center = FALSE,
            scale = FALSE)^2)
      }
      r.sq <- mss / (mss + rss)
      adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
      out <- list(R.squared = r.sq, adj.R.squared = adj.r.sq)
   } else {
      r.df <- summary(model)$df[2]
      int.df <- 1
      tss <- sum(w * (resp - center)^2)
      r.sq <- 1 - rss/tss
      adj.r.sq <- 1 - (1 - r.sq) * (n - int.df) / r.df
      out <- list(pseudo.R.squared = r.sq,
         adj.R.squared = adj.r.sq)
   }
   class(out) <- "Rsq"
   return(out)
}
