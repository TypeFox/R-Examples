simple.ef2 <- function(X, drift, sigma, h, h.x, h.xx,
                      guess, lower, upper){

 if(!is.ts(X))
  stop("Please provide a `ts' object")

 if(!is.expression(drift) | !is.expression(sigma))
  stop("Coefficients `drift' and `sigma' must be expressions") 	

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
 n.h <- length(h)
 if(!is.list(h) | (n.h!=n.vars-1) | 
   !all(unlist(lapply(h,mode))=="expression"))
  stop("`h' must be a list of expressions of length equal to the \
        dimension of the parameter space")
   
 new.env() -> e1
 
 self.hx <- FALSE
 if(missing(h.x)){
  message("h.x not provided, attempting symbolic derivation.\n")
  h.x <- vector(n.h, mode="list")
  for(i in 1:n.h)
   h.x[[i]] <- deriv(h[[i]],"x")
  self.hx <- TRUE
 }


 self.hxx <- FALSE
 if(missing(h.xx)){
  message("h.xx not provided, attempting symbolic derivation.\n")
  h.xx <- vector(n.h, mode="list")
  for(i in 1:n.h)
   h.xx[[i]] <- deriv(h[[i]],"x",hessian=TRUE)
  self.hxx <- TRUE
 }


 h1.x <- vector(n.h, mode="list")
 for(i in 1:n.h){
  h1.x[[i]] <- function(x,j) { 
   lx <- length(x)
   assign("x",x, e1)
   val <- eval(h.x[[j]],e1)
   if(self.hx)
    val <- as.numeric(attr(val,"gradient"))
   if(length(val) != lx)
    val <- rep(val, lx)
   val 
  }
 }
  

 h1.xx <- vector(n.h, mode="list")
 for(i in 1:n.h){
  h1.xx[[i]] <- function(x,j) { 
   lx <- length(x)
   assign("x",x, e1)
   val <- eval(h.xx[[j]],e1)
   if(self.hxx)
    val <- as.numeric(attr(val,"hessian")[1,1,1])
   if(length(val) != lx)
    val <- rep(val, lx)
   val 
  }
 }
 
 Dd <- deriv(drift, par.vars)
 Ds <- deriv(sigma, par.vars)

 D1 <- function(x){
  assign("x",x, e1)
  val <- as.numeric(eval(Dd,e1))
  if(is.na(d.has.x))
   val <- rep(val, length(x))
  val 
 }

 S1 <- function(x){
  assign("x",x, e1)
  val <- as.numeric(eval(Ds,e1))
  if(is.na(s.has.x))
   val <- rep(val, length(x))
  val 
 }

 Fn <- function(theta){
  for(i in 2:n.vars)
   assign(par.vars[i], theta[i-1],e1)	
  dd1 <-  D1(X)
  ss1 <- S1(X)
  val <- numeric(n.h)
  for(i in 1:n.h){
   hh1.x <-  h1.x[[i]](X,i)
   hh1.xx <- h1.xx[[i]](X,i)
   val[i] <- sum(dd1*hh1.x + 0.5*(ss1^2)*hh1.xx)
  }
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
 cat("\nOptimization constraints\n")
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
