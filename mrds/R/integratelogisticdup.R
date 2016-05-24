integratelogisticdup <- function(x1, x2, models, beta, lower=0, width, point){
  # numerical integral of product of logistic detection functions

  # computation speed-up when there is distance in the formula
  # but only if there are no interactions with distance
  if(sum(grepl("distance",names(beta)))==1){

    # set parameter to be zero for distance
    beta_distance <- beta[grepl("distance",names(beta))]
    beta[grepl("distance",names(beta))] <- 0

    # calculate the rest of the linear predictor
    x1 <- setcov(x1,models$g0model)%*%beta
    x2 <- setcov(x2,models$g0model)%*%beta

    # do some integration
    integrate(logisticdupbyx_fast, lower=lower, upper=width,
              subdivisions=10, rel.tol=0.01, abs.tol=0.01,
              x1=x1, x2=x2, models=models, beta=beta, point=point,
              beta_distance=beta_distance)$value
  }else{
  # Otherwise just go ahead and do the numerical integration
    integrate(logisticdupbyx, lower=lower, upper=width,
              subdivisions=10, rel.tol=0.01, abs.tol=0.01,
              x1=x1, x2=x2, models=models, beta=beta, point=point)$value
  }
}
