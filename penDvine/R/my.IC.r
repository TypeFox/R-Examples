my.IC <- function(penden.env,temp=FALSE) {
  values <- -eigen(get("Derv2.cal",penden.env))$values
  rho <- min(values[values>1e-09])
  d <- dim(get("Derv2.cal",penden.env))[1]
  if(temp) {
    mytrace <- sum(diag(solve(get("Derv2.pen.temp",penden.env),tol=1e-50)%*%(get("Derv2.cal.temp",penden.env))))
    assign("cAIC.temp",-get("log.like.temp",penden.env)+mytrace+(2*mytrace*(mytrace+1))/(get("n",penden.env)-mytrace-1),penden.env)
    assign("mytrace.temp",mytrace,penden.env)
  }
  else {
    mytrace <- sum(diag(solve(get("Derv2.pen",penden.env),tol=1e-50)%*%(get("Derv2.cal",penden.env))))
    assign("AIC",-get("log.like",penden.env)+mytrace,penden.env)
    assign("cAIC",-get("log.like",penden.env)+mytrace+(2*mytrace*(mytrace+1))/(get("n",penden.env)-mytrace-1),penden.env)
    assign("mytrace",mytrace,penden.env)
  }
}
