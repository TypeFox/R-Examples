linear.mart.ef <- function(X, drift, sigma, a1, a2, guess, 
  lower, upper, c.mean, c.var){

 if(!is.ts(X))
  stop("Please provide a `ts' object")

 N <- 1
 if(missing(a1))
  stop("`a1' is missing")
 if(!missing(a2)){
  if(length(a1) != length(a2))
   stop("`a1' and `a2' not of the same dimension")
  N <- 2 
 }

 if(!is.expression(a1))
  stop("`a1' must be a vector of expressions")

 if(N==2)
  if(!is.expression(a2))
  stop("`a2' must be a vector of expressions") 
 
 if(!is.expression(drift) | !is.expression(sigma))
  stop("Coefficients `drift' and `sigma' must be expressions") 	

 if(!is.expression(c.mean))
  stop("Conditional mean must be an expression")

 if(!is.expression(c.var))
  stop("Conditional variance must be an expression")
  
 d.vars <- all.vars(drift)
 s.vars <- all.vars(sigma)
 match("x", d.vars) -> d.has.x
 match("x", s.vars) -> s.has.x
 par.vars <- unique(c(d.vars, s.vars))
 n.vars <- length(par.vars)
 match("x", par.vars) -> idx
 if(is.na(idx))
  stop("One variable should be named `x'")
 par.vars <- par.vars[c(idx,(1:n.vars)[-idx])]
 n.pars <- n.vars - 1
 
# We check the list of expressions needed for the weights a(x,theta)
# each a_i(x,theta) must be a vector of expressions of the same 
# dimension of the parameter space
 a.vars <- all.vars(a1)
 if(N==2)
  a.vars <- c(a.vars, all.vars(a2))

 a.vars <- unique(a.vars)
 n.a.vars <- length(a.vars)
 match("x", a.vars) -> a.idx
 a.has.x <- !is.na(a.idx)
 a.idx <- as.integer(na.omit(a.idx))
 a.pars <- a.vars[(1:n.a.vars)[-a.idx]]
 n.a.pars <- length(a.pars)
 
 if(n.pars != length(a1))
  stop("weight functions `a' must have the same dimension of the\
  parameter space")

 new.env() -> e1
 # when missing c.mean and c.var need to be estimated
 # via MC methods. Easy.
 
 F.XT <- function(){
  val <- eval(c.mean,e1)
  if(length(val) != lx)
   val <- rep(val, lx)
  val  
 }

 PHI.XT <- function(){
  val <- eval(c.var,e1)
  if(length(val) != lx)
   val <- rep(val, lx)
  val  
 }
 
# vectorized version of a: A[i,j,] = a_{ij}() 
  A <- function() { 
   val <- array(0, c(N, n.pars, lx))
   for(k in 1:n.pars){  
	val[1,k,] <- eval(a1[k],e1)
	if(N==2)
	 val[2,k,] <- eval(a2[k],e1)
   }
   val 
  }
 

 n.obs <- length(X)
 Y.data <- X[2:n.obs]
 X.data <- X[1:(n.obs-1)]
  
  
 assign("x",X.data, e1)
 assign("y",Y.data, e1)
 lx <- length(X.data)

 Fn <- function(theta){
  for(i in 2:n.vars){
   assign(par.vars[i], theta[i-1],e1)	
   # assign(par.vars[i], theta[i-1], e1)
  }

  aa <-  A() # this contains the weights

  H1 <- Y.data - F.XT()
  H2 <- H1^2 -	PHI.XT()

  val <- 0
  for(i in 1:n.pars)
   val <- val + sum(aa[1,i,]*H1)
  if(N==2)
   for(i in 1:n.pars)
    val <- val + sum(aa[2,i,]*H2)
  return(val)
 }

 Gn <- function(theta)  sum(abs(Fn(theta))^2)
	

 if(missing(guess))
  start <- runif(n.vars-1)
 else
  start <- guess

 if(missing(lower))
  lower <- rep(-Inf, n.vars-1)

 if(missing(upper))
  upper <- rep(Inf, n.vars-1)

 st <- start
 names(st) <- par.vars[-1]
 cat("\nInitial values for the optimization algorithm ")
 if(missing(guess))
  cat("(random)\n")
 else
  cat("\n")
 print(st)
 cat("\nOptimization contraints\n")
 ct <- as.matrix(cbind(lower, upper))
 rownames(ct) <- par.vars[-1]
 colnames(ct) <- c("lower", "upper")
 print(ct)
  
 cat("\nRunning optimizer...\n")
 mn <- optim(start, Gn, method="L-BFGS-B", lower=lower, upper=upper)
 names(mn$par) <- par.vars[-1]

 estimate <- mn$par						
 fn <-  Fn(mn$par)

 return( list(estimate=estimate, Fn=fn) )
}