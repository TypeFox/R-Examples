
exact.midp=function(b,c,n,alpha){
		psidach=(b+c)/n
		tetadach=(b-c)/n
		binom.midp <- function(x, n, conf.level = 0.95) {
				  #Functions to find root of for the lower and higher bounds of the CI
				  #These are helper functions.
				  f.low <- function(pi,x,n) {
					1/2*dbinom(x,size=n,prob=pi) + pbinom(x,size=n,prob=pi,lower.tail=FALSE) - (1-conf.level)/2
				  }
				  f.up <- function(pi,x,n) {
					1/2*dbinom(x,size=n,prob=pi) + pbinom(x-1,size=n,prob=pi) - (1-conf.level)/2
				  }

				  #Function to calculate the midp confidence interval for one
				  #set of (x,n) values
				  one.ci <- function(x,n) {
				  #One takes pi_low = 0 when x=0 and pi_up=1 when x=n
					pi.low <- 0
					pi.up  <- 1

				  #Calculate CI by finding roots of the f funcs
					if (x!=0) {
					  pi.low <- uniroot(f.low,interval=c(0,x/n),x=x,n=n)$root
					} 
					if (x!=n) {
					  pi.up  <- uniroot(f.up,interval=c(x/n,1),x=x,n=n)$root
					}
					#Done
					return(c(pi.low,pi.up))
				  }

				  #Handle possible vector arguments
					#  if ((length(x) != length(n))) {
					m <- cbind(x = x, n = n)
					x <- m[, "x"]
					n <- m[, "n"]
					#  } 

				  #Run function
				  cis <- t(sapply(1:nrow(m), function(i) one.ci(x=x[i],n=n[i])))
				  #Wrap up result
				  res <- data.frame(method="midp",x=x,n=n,mean=x/n,lower=cis[,1],upper=cis[,2])

				  return(res)
				}
		midp=binom.midp(b,b+c,1-alpha) 
		exact2_low=(2*midp[[5]]-1)*psidach
		exact2_u=(2*midp[[6]]-1)*psidach
		cint=c(exact2_low,exact2_u)
		attr(cint, "conf.level") <- 1-alpha
		rval <- list(conf.int = cint, estimate = tetadach)
		class(rval) <- "htest"
		return(rval)
		}