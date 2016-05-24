library(mclogit)
options(error=recover)

mclogitP <- function(eta,s){
  expeta <- exp(eta)
  sum.expeta <- rowsum(expeta,s)
  expeta/sum.expeta[s]
}

N <- 10000
n <- 100

test.data <- data.frame(
  x = rnorm(N),
  f = gl(4,N/4),
  set = gl(N/5,5,N),
  altern0 = gl(5,1,N),
  nat = gl(15,N/15,N),
  occ = gl(10,1,N)
)

test.data <- within(test.data,{
  
    altern <- as.integer(interaction(altern0,nat))
    altern.occ <- as.integer(interaction(altern,occ))
    b1 <- rnorm(n=length(altern))
    b2 <- rnorm(n=length(altern.occ))
    ff <- 1+.2*(as.numeric(f)-1)
    eta <- x*ff + b1[altern] + b2[altern.occ]
    p <- mclogitP(eta,set)
    n <- unlist(tapply(p,set,function(p)rmultinom(n=1,size=n,prob=p)))
    rm(b1,b2)
})


test.mc0 <- mclogit(cbind(n,set)~x:f,data=test.data
              )


test.mc <- mclogit(cbind(n,set)~x:f,data=test.data,
              random=~1|altern/occ,
              #start.theta=c(1,1)
              maxit=100
              )

# By construction, the `true' coefficient values
# are 1, 1.2, 1.4, 1.6
coef(test.mc)

# The asymptotic covariance matrix of the coefficient estimates
vcov(test.mc)

print(test.mc)

summary(test.mc)

