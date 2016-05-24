#calculates the Dm-penalty matrix
D.m <- function(penden.env) {
 assign("Dm",kronecker(t(get("L",penden.env)),diag(1,get("p",penden.env)))%*%((kronecker(diag(1,(get("K",penden.env)-get("m",penden.env))),crossprod(get("Z",penden.env))/get("n",penden.env)))%*%kronecker(get("L",penden.env),diag(1,get("p",penden.env)))),penden.env)
}
