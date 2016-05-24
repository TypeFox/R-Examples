my.IC <- function(penden.env) {
  aa <- try(solve(get("Derv2.pen",penden.env)))
  eps <- 1e+08
  i <- 1
  while(class(aa)=="try-error") {
    aa <- try(solve(get("Derv2.pen",penden.env)+diag(i*eps,dim(get("Derv2.pen",penden.env))[1])))
    i <-  i+1
  }
  if(get("base",penden.env)=="B-spline") mytrace <- sum(diag(aa%*%(get("Derv2.cal",penden.env))))#-get("linear.con",penden.env)
  if(get("base",penden.env)=="Bernstein") mytrace <- get("DD",penden.env)#-get("linear.con",penden.env)

  #if(get("base",penden.env)=="B-spline") {
    #diff <- get("DD",penden.env)-dim(get("AA.help",penden.env))[2]
    #diff3 <- dim(get("AA.help",penden.env))[2]+1
    #C <- t(cbind(get("AA.help",penden.env),rbind(diag(val,diff),cbind(diag(val,diff2),matrix(rep(0,diff2*diff3),ncol=(diff3))))))
    #C <- t(cbind(get("AA.help",penden.env),diag(1,get("DD",penden.env))[,c(diff3:get("DD",penden.env))]))

    #C <- t(cbind(get("AA.help",penden.env),diag(c(1/(get("ck.val",penden.env))))))[seq(1,get("DD",penden.env)),seq(1,get("DD",penden.env))]
   
   # H1 <- (t(my.positive.definite.solve(C))%*%get("Derv2.pen",penden.env)%*%my.positive.definite.solve(C))[seq(diff3,get("DD",penden.env)),seq(diff3,get("DD",penden.env))]
   # H2 <- (t(my.positive.definite.solve(C))%*%get("Derv2.cal",penden.env)%*%my.positive.definite.solve(C))[seq(diff3,get("DD",penden.env)),seq(diff3,get("DD",penden.env))]
   # mytrace <- sum(diag(solve(H1)%*%H2))
  #}

  #return(list(myAIC=(-2*get("log.like",penden.env)+2*mytrace),mytrace=mytrace,myBIC=(2*get("log.like",penden.env)+mytrace*log(get("n",penden.env)))))
  assign("AIC",-2*get("log.like",penden.env)+2*mytrace,penden.env)
  assign("cAIC",get("AIC",penden.env)+(2*mytrace*(mytrace+1))/(get("n",penden.env)-mytrace-1),penden.env)
  assign("BIC",-2*get("log.like",penden.env)+mytrace*log(get("n",penden.env)),penden.env)
  assign("mytrace",mytrace,penden.env)
  # if(get("base",penden.env)=="B-spline")assign("mytrace2",mytrace2,penden.env)
}
