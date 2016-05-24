nnls2 <- function (cvec, beta.hat.null, x.tilde, y, family, lambda.nng, nvars)
{
     eta.new <- drop (beta.hat.null + x.tilde %*% cvec)
     mu.new <- family$linkinv (eta.new)      # fitted values
     d.new <- family$mu.eta (eta.new)        # derivative of response function
     v.new <- family$variance (mu.new)       # variance function of the response
     weights <- d.new / sqrt (v.new)  # decomposed elements (^0.5) of weight matrix W, see GLM notation
     x.star <- weights * x.tilde  
     y.tilde.star <- weights * (eta.new  + (y - mu.new) / d.new)    
     
     nnls.y <- c (y.tilde.star, rep (0, nvars))
     nnls.x <- rbind (x.star, matrix (sqrt (lambda.nng), nrow = nvars, ncol = nvars))
  
     return (sum ((drop (nnls.y - nnls.x %*% cvec))^2))

}
