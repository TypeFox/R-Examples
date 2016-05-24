simple.ef <- function(X, f, guess, lower, upper){

 if(!is.ts(X))
  stop("Please provide a `ts' object")

 n.f <- length(f)
 f.vars <- NULL
 for(i in 1:n.f)
  f.vars <- c(f.vars,all.vars(f[[i]]))
 f.vars <- unique(f.vars)
 n.vars <- length(f.vars)
 match("x", f.vars) -> idx.x
 match("y", f.vars) -> idx.y
 has.x <- !is.na(idx.x)
 has.y <- !is.na(idx.y)
 if(!has.x){
  if(!has.y)
   stop("Variables `x' and `y' both missing")
  else
   stop("Variable `x' missing")   
 } 
 idx <- c(idx.x, idx.y) 
 idx <- as.integer(na.omit(idx))
 f.pars <- f.vars[(1:n.vars)[-idx]]
 n.pars <- length(f.pars)
 f.vars <- f.vars[idx]
 
 if(!is.list(f) | (n.f!=n.pars) | 
   !all(unlist(lapply(f,mode))=="expression"))
  stop("`f' must be a list of expressions of length equal to the \
        dimension of the parameter space")
   
 new.env() -> e1

 f1 <- vector(n.f, mode="list")
 for(i in 1:n.f){
  f1[[i]] <- function(j) { 
   val <- eval(f[[j]],e1)
   if(length(val) != lx)
    val <- rep(val, lx)
   val 
  }
 }

 X.data <- NA
 Y.data <- NA 
 n.obs <- length(X)
 if(has.y)
  Y.data <- X[2:n.obs]
 if(has.x){
  if(has.y)
   X.data <- X[1:(n.obs-1)]
  else
   X.data <- X[1:n.obs]
 }
 assign("x",X.data, e1)
 assign("y",Y.data, e1)
 lx <- length(X.data)
 
 Fn <- function(theta){
  for(i in 1:n.pars)
   assign(f.pars[i], theta[i], e1)	
  val <- numeric(n.f)
  for(i in 1:n.f){
   ff1 <-  f1[[i]](i)
   val[i] <- sum(ff1)
  }
  return(val)
 }

 Gn <- function(theta)  sum(abs(Fn(theta))^2)
	

 if(missing(guess))
  start <- runif(n.pars)
 else
  start <- guess

 if(missing(lower))
  lower <- rep(-Inf, n.pars)

 if(missing(upper))
  upper <- rep(Inf, n.pars)

 st <- start
 names(st) <- f.pars
 cat("\nInitial values for the optimization algorithm ")
 if(missing(guess))
  cat("(random)\n")
 else
  cat("\n")
 print(st)
 cat("\nOptimization constraints\n")
 ct <- as.matrix(cbind(lower, upper))
 rownames(ct) <- f.pars
 colnames(ct) <- c("lower", "upper")
 print(ct)
  
 cat("\nRunning optimizer...\n")
 mn <- optim(start, Gn, method="L-BFGS-B", lower=lower, upper=upper)
 names(mn$par) <- f.pars

 estimate <- mn$par						
 fn <-  Fn(mn$par)

 return( list(estimate=estimate, Fn=fn) )
}
