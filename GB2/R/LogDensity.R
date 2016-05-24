
logf.gb2 <- function(x, shape1, scale, shape2, shape3){
y <- (x/scale)^shape1
logf <- log(shape1/scale) - lbeta(shape2,shape3) + (shape1*shape2 - 1)*log(x/scale) - (shape2+shape3)*log(1+y)
return(logf)
}

dlogf.gb2 <- function(xi, shape1, scale, shape2, shape3){
y <- (xi/scale)^shape1
logy <- log(y)
dlogf.dshape1 <- (1/shape1)*(1 + shape2*logy - (shape2+shape3)*y*logy/(1+y))
dlogf.dscale <- -(shape1/scale)*shape2 + (shape1/scale)*(shape2+shape3)*y/(1+y)
dlogf.dshape2 <- digamma(shape2+shape3) - digamma(shape2) + logy - log(1+y)
dlogf.dshape3 <- digamma(shape2+shape3) - digamma(shape3) - log(1+y)
dlogf <- c(dlogf.dshape1, dlogf.dscale, dlogf.dshape2, dlogf.dshape3)
return(dlogf)
}

d2logf.gb2 <- function(xi, shape1, scale, shape2, shape3){
y <- (xi/scale)^shape1
logy <- log(y)
h1 <- y/(y+1)
h2 <- logy/(y+1)
h3 <- h1*h2
h4 <- h3*logy
trip <- trigamma(shape2)
triq <- trigamma(shape3)
trippq <- trigamma(shape2+shape3)
D <- matrix(rep(0,16), ncol=4)
D[1,1]  <- -1/(shape1^2) - ((shape2+shape3)*h4)/(shape1^2)
D[1,2]  <- -shape2/scale + (shape2+shape3)*(h1 + h3)/scale
D[1,3]  <- h2/shape1
#D[1,4] <- h2/shape1 - logy/shape1
D[1,4]  <- (-logy/shape1)*h1
D[2,1]  <- D[1,2]
D[2,2]  <- (shape1*shape2)/(scale^2) - shape1*(shape2+shape3)*(h1 + shape1*h1/(y+1))/(scale^2)
#D[2,3] <- (shape1/scale)*h1 - shape1/scale
D[2,3]  <- -shape1/(scale*(y+1))
D[2,4] <- (shape1/scale)*h1
D[3,1] <- D[1,3]
D[3,2] <- D[2,3]
D[3,3] <- trippq - trip
D[3,4] <- trippq
D[4,1] <- D[1,4]
D[4,2] <- D[2,4]
D[4,3] <- D[3,4]
D[4,4] <- trippq - triq
return(D)
}
