ellip <- function(type = c("chi","t2"),x, Xmv, S, phase = 1, alpha = 0.01,method = "sw",colm,  ...){

type <- match.arg(type)

p <- ncol(x) # quality characteristics
m <- nrow(x) # number of samples or observations
if (class(x) == "matrix" || class(x) == "data.frame") (x <- array(data.matrix(x),c(m,p,1)))
n <- dim(x)[3] # observations or sample size 
 
if(!missing(Xmv))(phase <- 2)

x.jk <- matrix(0,m,p)

x.jk <- apply(x,1:2,mean)

if(missing(Xmv))(Xmv <- colMeans(x.jk)) 
if(missing(S))(S <- covariance(x,method = method))
if(missing(colm))(colm <- nrow(x))


Ue <- eigen(S)$vectors  
DDe <- eigen(S)$values



if (type == "chi") { di <- sqrt(qchisq(1 - alpha, p)/n)}
else {
ifelse(n == 1, ifelse(phase == 1, 
 di <- sqrt(((colm - 1) ^ 2) / colm * qbeta(1 - alpha,p / 2,(((2 * (colm - 1) ^ 2) / (3 * colm - 4) - p - 1) / 2))),
 di<-sqrt(((p * (colm + 1) * (colm - 1)) / ((colm ^ 2) - colm * p)) * qf(1 - alpha,p,colm - p))),
 ifelse(phase == 1, 
 di <- sqrt((p * (colm - 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1)/n),
 di <- sqrt((p * (colm + 1) * (n - 1)) / (colm * n - colm - p + 1) * qf(1 - alpha,p,colm * n - colm - p + 1)/n))
)

}

angle <- seq(0, 2 * pi, length.out = 200)
ch1 <- cbind(di * cos(angle),di * sin(angle))
tt <- t(Xmv - ((Ue %*% diag(sqrt(DDe))) %*% t(ch1)))

     return(tt)
}
