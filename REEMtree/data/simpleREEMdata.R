# This creates a simple example that can be used for the R packages


set.seed(39)

T <- 12
K <- 50
timeoffset <- 1
dummyoffset <- 2

alpha <- rnorm(K, sd=2)
eps <- rnorm(T*K)

data <- matrix(nrow=T*K,ncol=5)
# Covariate is time
data[,2] <- 1:T
# Other covariate is randomly generated
data[,5] <- rnorm(T*K, sd=0.5)
data[,3] <- c(rep(0,T*K/2),rep(1,T*K/2))
for(i in 1:K){
    data[((i-1)*T+1):(i*T),1] <- alpha[i] + 
	dummyoffset*(i>K/2)+
	timeoffset*(data[((i-1)*T+1):(i*T),2]>T/2) *(i>K/2)+
	(data[((i-1)*T+1):(i*T),5]<0.25)*(i<=K/2)+
	eps[((i-1)*T+1):(i*T)]
    data[((i-1)*T+1):(i*T),4] <- i
}

simpleREEMdata <- data.frame(data)
names(simpleREEMdata) <- c("Y","t","D", "ID", "X")

# Removing unnecessary variables
rm(T, K, timeoffset, dummyoffset, alpha, eps, i, data)

