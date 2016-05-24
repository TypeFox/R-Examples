`ss.aipe.sem.path.sensitiv` <-
function(model, est.Sigma, true.Sigma=est.Sigma, which.path,
desired.width, N=NULL, conf.level=0.95, assurance=NULL, G=100, ...){

  if(!requireNamespace("MASS", quietly = TRUE)) stop("The package 'MASS' is needed; please install the package and try again.")
  if(!requireNamespace("sem", quietly = TRUE)) stop("The package 'sem' is needed; please install the package and try again.")
  
result.plan <- ss.aipe.sem.path(model=model, Sigma=est.Sigma, desired.width=desired.width, 
which.path=which.path, conf.level=conf.level, assurance=assurance, ...)

obs.vars<- result.plan$obs.vars
J <- length(result.plan$parameter)
alpha <- 1-conf.level
p <- dim(est.Sigma)[1]
if(is.null(N)) N<- result.plan$sample.size
j <- result.plan$path.index

Data<- matrix(NA, N,p)
S<- matrix(NA, p,p)
theta.hat.j <- rep(NA, G)
#theta.hat <- array(NA, c(J,1,G))
SE.theta.hat.j <- rep(NA, G)
cov.theta.hat <- matrix(NA, J,J)

g<- 0
while(g<G){
	gc()
	cat("successful.iteration = ", g, "\n")
   Data <- MASS::mvrnorm(n = N, mu=rep(0,p), Sigma=true.Sigma)
   S <- var(Data)  
   colnames(S) <- rownames(S)<- obs.vars
   S.fit <- try(sem::sem(model, S, N), FALSE) 
   if(inherits(S.fit, "try-error")|| S.fit$convergence>2) g<- g
   else{
   		g<- g+1
   		cov.theta.hat <- S.fit$cov
      SE.theta.hat.j[g] <- ifelse (any(diag(cov.theta.hat)< 0), NA, sqrt(cov.theta.hat[j,j]))
      theta.hat.j[g] <- S.fit$coeff[j]
   		} 
	}#end of while(suc.rep<G)
	
    
CI.upper <- theta.hat.j+ qnorm(1-alpha/2)*SE.theta.hat.j
CI.lower <- theta.hat.j- qnorm(1-alpha/2)*SE.theta.hat.j
w <- CI.upper - CI.lower

true.fit <- sem::sem(model, true.Sigma, 100000)
theta.j <- true.fit$coeff[j]

CI.upper <- na.omit(CI.upper)
CI.lower <- na.omit(CI.lower)
w <- na.omit(w)
G<- length(w)

percent.narrower <- sum(w<=desired.width)/G
alpha.emp.upper <- sum(theta.j>CI.upper)/G
alpha.emp.lower <- sum(theta.j<CI.lower)/G

result<- list()
result$w <- w
result$sample.size <- N
result$path.of.interest <- which.path
result$desired.width <- desired.width
result$mean.width <- mean(w)
result$median.width <- median(w)
result$quantile.width <- quantile(w, c(.99,.95,.90,.85,.80, .75, .70, .60))
result$width.less.than.desired <- percent.narrower
result$Type.I.err.upper <- alpha.emp.upper
result$Type.I.err.lower <- alpha.emp.lower
result$Type.I.err <- alpha.emp.upper+alpha.emp.lower
result$conf.level <- conf.level
result$rep <- G
return(result)  
}# end of function()

