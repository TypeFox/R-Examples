ran.data <-
function (d, t, mu=1, lambda=5, A=15)
{
many.data  <- d[1]
many.time  <- t[1]
mu         <- mu[1]
lambda     <- lambda[1]
A          <- A[1]

# /// check input
if (is.numeric(many.data)==FALSE) stop("Need numeric value for: d")
if (many.data<30) stop("Need number of datasets d > 30")
if (is.numeric(many.time)==FALSE) stop("Need numeric value for: t")
if (is.numeric(mu)==FALSE) stop("Need numeric value for: mu")
if (is.numeric(lambda)==FALSE) stop("Need numeric value for: lambda")
if (is.numeric(A)==FALSE) stop("Need numeric value for: A")

# /// set up time + noise
time  <- t(array(rep(c(1:many.time),many.data),c(many.time,many.data)))+ runif(many.data*many.time)

# /// set up growth data array
data  <- data.frame(array(0,c(many.data,(many.time+3))))

# /// number of datasets per test
nT1   <- round(d/3) 
nT2   <- round(d/3) 
nT3   <- d - nT1 - nT2
test1 <- rep("Test I"  , nT1)
test2 <- rep("Test II" , nT2)
test3 <- rep("Test III", nT3)

# /// concentrations + noise
conc1 <- seq(0.1, 1, length.out=nT1) + runif(nT1, -0.05, 0.05)
conc2 <- seq(0.1, 1, length.out=nT2) + runif(nT2, -0.05, 0.05)
conc3 <- seq(0.1, 1, length.out=nT3) + runif(nT3, -0.05, 0.05)

data[,1] <- (c(test1, test2, test3))
data[,3] <- (c(conc1, conc2, conc3))

# /// set parameter changes
mu.loop     <- mu - mu * data[,3]^2/( (mu/2)^2 + data[,3]^2 )
A.loop      <- A - 0.75 * A * data[,3]
lambda.loop <- lambda + 0.5 * lambda * data[,3]/(0.5+data[,3])

# /// create growth data
for (i in 1:many.data){
data[i,4:(many.time+3)] <- gompertz(time[i,], A.loop[i], mu.loop[i], lambda.loop[i]) + 0.02*A*rnorm(many.time)
}

# /// add information columns
data     <- data.frame(data)

data[,2] <- (rep(c("A", "B", "C", "D", "E", "F", "G"), many.data ))[1:many.data]

ran.data <- list(data=data, time=time)

}

