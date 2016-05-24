#Calcualtes the AIC-Criteria for the estimated density.
my.AIC <- function(penden.env,lambda0,opt.Likelihood=NULL) {
  if(!is.null(opt.Likelihood)) {val1 <- -opt.Likelihood}
  if(is.null(opt.Likelihood)) {val1 <- -pen.log.like(penden.env,lambda0=0)}
  df <- my.positive.definite.solve(get("Derv2.pen",penden.env))%*%get("Derv2.cal",penden.env)
  mytrace <- sum(diag(df))
  return(list(myAIC=(val1+mytrace),mytrace=mytrace))
}
