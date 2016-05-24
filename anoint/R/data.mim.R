# SIMULATE LOGISTIC PIM DATASET

data.anoint <- function(

				alpha,
				beta,
				gamma,
				mean,
				vcov,
				n = 100,
				event = .8,
				type = c("binomial","survival")
)
{

data.logit.anoint <- function(

				alpha,
				beta,
				gamma,
				mean,
				vcov,
				n = 100	
)
{
	
	trt <- rep(c(1,0),each=n)
	X0 <- rmnorm(n=n,mean,vcov)
	X1 <- rmnorm(n=n,mean,vcov)
	
	expit <- function(x) exp(x)/(1+exp(x))

	pi0 <- expit(alpha[1]+X0%*%beta)
	pi1 <- expit(alpha[2]+X1%*%diag(gamma)%*%beta)
	
	y <- sapply(c(pi1,pi0),function(p)rbinom(n=1,size=1,prob=p))
	
	df <- as.data.frame(rbind(X1,X0))
	df$y <- y
	df$trt <- trt
	
 df
}


data.cox.anoint <- function(

				alpha,
				beta,
				gamma,
				mean,
				vcov,
				n = 100,
				event = .8
)
{
	
	trt <- rep(c(1,0),each=n)
	X0 <- rmnorm(n=n,mean,vcov)
	X1 <- rmnorm(n=n,mean,vcov)

	#SPECIFY RATE OF EXPONENTIAL
	r0 <- exp(alpha[1]+X0%*%beta) 
	r1 <- exp(alpha[2]+X1%*%diag(gamma)%*%beta)
	
	y <- sapply(c(r1,r0),function(p)rexp(n=1,rate=p)) # TIMES
	
	df <- as.data.frame(rbind(X1,X0))
	df$y <- y
	
	# DETERMINE OBSERVED EVENTS
	df$trt <- trt
	df$event <- 1
	
	# FOR NON-INFORMATIVE RIGHT-CENSORED
	censored <- sample(nrow(df),ceiling((1-event)*nrow(df)))
	# TAKE RANDOM EXIT FROM OBSERVATION
	exit <- sapply(df$y[censored],function(t)runif(n=1,0,max=t))
	df$y[censored] <- exit
	df$event[censored] <- 0
	
 df
}

	if(type=="survival")
		data.cox.anoint(alpha,beta,gamma,mean,vcov,n,event)	
	else
		data.logit.anoint(alpha,beta,gamma,mean,vcov,n)

}
