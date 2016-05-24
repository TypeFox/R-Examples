biascorr <-
function(b, b.se, p.alpha, T){

  c_alpha <- abs(qt((1-(1-p.alpha))/2, T))
  bt <- b/b.se

  dr <- (dnorm(c_alpha-bt)-dnorm(-c_alpha-bt))/ (1-pnorm(c_alpha-bt) + pnorm(-c_alpha-bt))
  dtbar <- bt - dr
  drbar <- (dnorm(c_alpha-dtbar)-dnorm(-c_alpha-dtbar))/ (1-pnorm(c_alpha-dtbar) + pnorm(-c_alpha-dtbar))

  #only correct if significant

  b_1step <- b
  b_1step[abs(bt)>c_alpha] <- b*(1-(dr/bt))[abs(bt)>c_alpha]

  b_2step <- b
  b_2step[abs(bt)>c_alpha] <- b*(1-(drbar/bt))[abs(bt)>c_alpha]

  b_corr <- cbind(b, b_1step, b_2step)
  colnames(b_corr) <- c("beta", "beta.1step", "beta.2step")

  return(b_corr)
}
