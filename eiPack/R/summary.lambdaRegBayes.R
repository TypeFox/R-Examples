summary.lambdaRegBayes <- function(object, ...){
  if(!is.mcmc(object$lambda.out)){
    coef <- apply(object$lambda.out, c(1,2), mean)
    sd <- apply(object$lambda.out, c(1,2), sd)
    quants <- matrix(apply(object$lambda.out, c(1,2), quantile, probs=c(.025,
                                                           0.05, 0.25,
                                                           0.5, 0.75,
                                                           0.95,
                                                           0.975)),
                   nrow(coef)*ncol(coef), 7, byrow=T)
    nidx <- apply(expand.grid(dimnames(coef)), 1, paste, collapse =
                  ".")
    tab <- cbind(c(coef), c(sd), quants)
    rownames(tab) <- nidx
  }
  if(is.mcmc(object$lambda.out)){
    coef <- apply(object$lambda.out, 2, mean)
    sd <- apply(object$lambda.out, 2, sd)
    quants <- t(apply(object$lambda.out, 2, quantile, probs=c(.025,
                                                           0.05, 0.25,
                                                           0.5, 0.75,
                                                           0.95,
                                                           0.975)))
    tab <- cbind(c(coef), c(sd), quants)
  }

  colnames(tab) <- c("Mean", "Std. Dev.", "2.5%", "5%", "25%",
                     "50%", "75%", "95%", "97.5%")
  tab
}
