
dgb2 <- function(x,shape1,scale,shape2,shape3){
	y <- (x/scale)^shape1
	dy_dx <- (shape1/scale)*(x/scale)^(shape1-1)
	z <- y/(1+y)
	z[z==1] <- 1-.Machine$double.eps
 	dz_dy <- (1+y)^(-2)
	dens <- dbeta(z, shape2, shape3) * dz_dy * dy_dx
    v <- (x==Inf)          
    dens[v] <- 0           
    return(dens)           
}

pgb2 <- function(x,shape1,scale,shape2,shape3){
	y <- (x/scale)^shape1
	z <- y/(1+y)
	prob <- pbeta(z, shape2, shape3)
    v <- (x==Inf)         
    prob[v] <- 1  
    return(prob)
}

qgb2 <- function (prob, shape1, scale, shape2, shape3) 
{
    pr <- sort(prob)
    ord <- order(prob)
    z1 <- qbeta( pr[pr<=0.5], shape2, shape3)
    z2 <- qbeta(1-pr[pr>0.5], shape3, shape2)
    y <- c( z1/(1 - z1), (1-z2)/z2 )
    quant <- y[ord]
    return(scale * quant^(1/shape1))
}

rgb2 <- function(n,shape1,scale,shape2,shape3){
	z <- rbeta(n,shape2,shape3)
	y <- z/(1-z)
	return(scale*y^(1/shape1))
}

