if(!require("GNE"))stop("this test requires package GNE.")

#-------------------------------------------------------------------------------
# (2) Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------


#constants
myarg <- list(d= 20, lambda= 4, rho= 1)

dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	res <- -arg$rho * x[i]
	if(i == j)
		res <- res + arg$d - arg$lambda - arg$rho*(x[1]+x[2])
	-res
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
	arg$rho * (i == j) + arg$rho * (j == k)	


dimlam <- c(1, 1)
#constraint function g_i(x)
g <- function(x, i)
	-x[i]
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	-1*(i == j)
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k)
	0

#true value is ((d-lambda)/3*rho, (d-lambda)/3*rho, 0, 0) 


getNE <- function(myarg, control=list(itermax=100, trace=0))
{
	z0 <- rep(0, sum(dimx)+sum(dimlam))
	GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=control$trace, maxit=control$itermax))$par[1:length(dimx)]
}
getNE4sobol <- function(x, control=list(itermax=100, trace=0))
{
#	res <- vector("numeric", NROW(x))
	print(dim(x))
	
	#for(i in 1:NROW(x))
	res <- sapply(1:NROW(x), function(i) 
	{		
		myarg <- list(d=x[i,1], lambda=x[i,2], rho=x[i,3])
		z0 <- rep(0, sum(dimx)+sum(dimlam))
		GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
			constr=g, NULL, grconstr=grg, NULL, heconstr=heg, NULL, 
			compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
			control=list(trace=control$trace, 
			maxit=control$itermax))$par[1]
	}
	)
	res
}

checkNE <- function(myarg)
{
	rep((myarg$d-myarg$lambda)/(3*myarg$rho), 2)
}

getNE(myarg)
checkNE(myarg)

sensiNE <- function(x, myarg, echo=FALSE)
{
#	x <- match.arg(x, c("d", "lambda", "rho"))
	for(v in x)
	{
		myarg[[v]] <- myarg[[v]] * (1+rnorm(1, 0, .1))
	}
	if(echo)
		print(unlist(myarg))
	
	c(unlist(myarg), ne=getNE(myarg))
}
rsensiNE <- function(nbsimu, x, myarg)
	t(replicate(nbsimu, sensiNE(x, myarg)))


sensid <- rsensiNE(10, "d", myarg)
plot(sensid[,"d"], sensid[,"ne1"])

sensirho <- rsensiNE(10, "rho", myarg)
plot(sensid[,"rho"], sensid[,"ne1"])

sensiall <- rsensiNE(10, c("d", "lambda", "rho"), myarg)

par(mfrow=c(1,3))
plot(sensiall[,"d"], sensiall[,"ne1"])
plot(sensiall[,"rho"], sensiall[,"ne1"])
plot(sensiall[,"lambda"], sensiall[,"ne1"])


y <- sensiall[,"ne1"]
shat <- sd(y)
x1 <- sensiall[,"d"]
shat1 <- sd(x1)
x2 <- sensiall[,"lambda"]
shat2 <- sd(x2)
x3 <- sensiall[,"rho"]
shat3 <- sd(x3)

# standardized regression coefficients
abs( coef(lm(y~x1))["x1"] * shat1 / shat )
abs( coef(lm(y~x2))["x2"] * shat1 / shat )
abs( coef(lm(y~x3))["x3"] * shat1 / shat )

coef(lm(y~x1+x2+x3))["x1"] * shat1 / shat

# data <- data.frame(Y=y, data.frame(x1, x2, x3))
# i <- 1:nrow(data)

    # d <- data[i, ]
    # lm.Y <- lm(formula(paste(colnames(d)[1], "~", paste(colnames(d)[-1], 
        # collapse = "+"))), data = d)
    # coefficients(lm.Y)[-1] * sapply(d[-1], sd)/sapply(d[1], sd)



#R squared
summary(lm(y~x1))$r.squared
summary(lm(y~x2))$r.squared
summary(lm(y~x3))$r.squared

summary(lm(y~x3))$r.squared
summary(lm(y~x3+x2))$r.squared
summary(lm(y~x3+x2+x1))$r.squared

#partial correlation
yhat1 <- fitted(lm(y~x2+x3))
x1hat <- fitted(lm(x1~x2+x3))

cor(x1-x1hat, y-yhat1)
cor(x1-x1hat, y-yhat1, m="kendall")
cor(x1-x1hat, y-yhat1, m="spearman")

yhat2 <- fitted(lm(y~x1+x3))
x2hat <- fitted(lm(x2~x1+x3))

cor(x2-x2hat, y-yhat2)
cor(x2-x2hat, y-yhat2, m="kendall")
cor(x2-x2hat, y-yhat2, m="spearman")

yhat3 <- fitted(lm(y~x1+x2))
x3hat <- fitted(lm(x3~x1+x2))

cor(x3-x3hat, y-yhat3)
cor(x3-x3hat, y-yhat3, m="kendall")
cor(x3-x3hat, y-yhat3, m="spearman")

par(mfrow=c(1,3))
plot(x1-x1hat, y-yhat1)
plot(x2-x2hat, y-yhat2)
plot(x3-x3hat, y-yhat3)


# packages
library(sensitivity)

#partial correlation coef
respcc <- pcc(data.frame(x1, x2, x3), y)
respcc$PCCC

#partial rank correlation coef
resprc <- pcc(data.frame(x1, x2, x3), y, rank=TRUE)
resprc$PRCC

ressrc <- src(data.frame(x1, x2, x3), y)
ressrc

ressrrc <- src(data.frame(x1, x2, x3), y, rank=TRUE)
ressrrc

n <- 10
X1 <- data.frame(x1 = runif(n, 10, 30), x2 = runif(n, 2, 6), x3 = runif(n, 1/2, 3))
X2 <- data.frame(x1 = runif(n, 10, 30), x2 = runif(n, 2, 6), x3 = runif(n, 1/2, 3))


system.time(res1 <- sobol(getNE4sobol, X1, X2, order=2))

plot(res1, ylim=c(-.25, 1), main="Sobol indices")

res12002 <- sobol2002(getNE4sobol, X1, X2)
res12002
res12007 <- sobol2007(getNE4sobol, X1, X2)
res12007




