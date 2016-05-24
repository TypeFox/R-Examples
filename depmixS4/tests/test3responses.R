# 
# test response models
# 

require(depmixS4)

# binomial response model
x <- rnorm(1000)
# library(boot)
invlogit <- function(lp) {exp(lp)/(1+exp(lp))}
p <- invlogit(x)
ss <- rbinom(1000,1,p)
mod <- GLMresponse(cbind(ss,1-ss)~x,family=binomial())
fit(mod)
glm(cbind(ss,1-ss)~x, family=binomial)

# gamma response model
x=runif(1000,1,5)
res <- rgamma(1000,x)
## note that gamma needs proper starting values which are not
## provided by depmixS4 (even with them, this may produce warnings)
mod <- GLMresponse(res~x,family=Gamma(),pst=c(0.8,1/0.8))
fit(mod)
glm(res~x,family=Gamma)

# multinomial response model
x <- sample(0:1,1000,rep=TRUE)
mod <- GLMresponse(sample(1:3,1000,rep=TRUE)~x,family=multinomial(),pstart=c(0.33,0.33,0.33,0,0,1))
mod@y <- simulate(mod)
fit(mod)
colSums(mod@y[which(x==0),])/length(which(x==0))
colSums(mod@y[which(x==1),])/length(which(x==1))

# multivariate normal response model
mn <- c(1,2,3)
sig <- matrix(c(1,.5,0,.5,1,0,0,0,2),3,3)
y <- mvrnorm(1000,mn,sig)
mod <- MVNresponse(y~1)
fit(mod)
colMeans(y)
var(y)

# normal (gaussian) response model
y <- rnorm(1000)
mod <- GLMresponse(y~1)
fm <- fit(mod)
cat("Test gaussian fit: ", all.equal(getpars(fm),c(mean(y),sd(y)),check=FALSE))


# poisson response model
x <- abs(rnorm(1000,2))
res <- rpois(1000,x)
mod <- GLMresponse(res~x,family=poisson())
fit(mod)
glm(res~x, family=poisson)



