if(!require("GNE"))stop("this test requires package GNE.")

#-------------------------------------------------------------------------------
# (2) generalized Duopoly game of Krawczyk and Stanislav Uryasev (2000)
#-------------------------------------------------------------------------------

obj <- function(x, i, d, lam, rho)
	(d - lam - rho*sum(x))*x[i]
	
NIF <- function(x, y, par)	
	obj(x, 1, par["d1"], par["lam1"], par["rho1"]) - obj(c(y[1], x[2]), 1, par["d1"], par["lam1"], par["rho1"]) + obj(x, 2, par["d2"], par["lam2"], par["rho2"]) - obj(c(x[1], y[2]), 2, par["d2"], par["lam2"], par["rho2"])
	
grobj <- function(x, i, j, d, lam, rho)
{
	if(i == j)
		d - lam - rho*sum(x) - rho*x[i]
	else
		-rho*x[j]
}	

#
getNE <- function(param)
{
	x1star <- 2*(param["d1"] - param["lam1"])/3/param["rho1"] - (param["d2"] - param["lam2"])/3/param["rho2"]
	x2star <- 2*(param["d2"] - param["lam2"])/3/param["rho2"] - (param["d1"] - param["lam1"])/3/param["rho1"]
	c(x1star, x2star)
}



getNE(c(d1=20, d2=20, lam1=4, lam2=4, rho1=1, rho2=1)) #(16/3, 16/3)

getNE(c(d1=20, d2=22, lam1=4, lam2=5, rho1=1, rho2=1.5))


genpar <- function()
{
	d <- runif(2, 15, 25)
	names(d) <- paste("d", 1:2, sep="")
	lam <- runif(2, 3, 5)
	names(lam) <- paste("lam", 1:2, sep="")
	rho <- runif(2, 1, 3)
	names(rho) <- paste("rho", 1:2, sep="")
	c(d, lam, rho)
}

#--------------------------------
#one-player sensi analysis

getsensiobj <- function(param, player=1)
{
	sapply(1:NROW(param), function(i)
		{
			par <- param[i, ]
			#print(par)
			xstar <- getNE(par)
			as.numeric( obj(xstar, player, par[paste("d", player, sep="")],
				par[paste("lam", player, sep="")],
				par[paste("rho", player, sep="")]) )
		}
	)
	
}

library(sensitivity)

n <- 100
X1 <- t(replicate(n, genpar()))
X2 <- t(replicate(n, genpar()))
Y1 <- getsensiobj(X1)

d1 <- data.frame(y=Y1, X1)

# plot(y~., data=d1)
# plot(d1)

scr1 <- src(X1, Y1)
pcc1 <- pcc(data.frame(X1), Y1)
pcc1$PCC

r1 <- sobol(getsensiobj, X1, X2, order=1)


playcomp <- cbind(scr1$SRC, pcc1$PCC, r1$S)
colnames(playcomp) <- paste("Play. Comp.", c("SCR", "PCC", "Sobol"), sep="-")


getsensiNE <- function(param, player=1)
{
	sapply(1:NROW(param), function(i)
		{
			par <- param[i, ]
			#print(par)
			xstar <- getNE(par)
			xstar[player]
		}
	)
	
}

X1 <- t(replicate(n, genpar()))
X2 <- t(replicate(n, genpar()))
Y1 <- getsensiNE(X1)

scr2 <- src(X1, Y1)
pcc2 <- pcc(data.frame(X1), Y1)
pcc2$PCC

r2 <- sobol(getsensiNE, X1, X2, order=1)

playobj <- cbind(scr2$SRC, pcc2$PCC, r2$S)
colnames(playobj) <- paste("Play. Obj.", c("SCR", "PCC", "Sobol"), sep="-")


#--------------------------------
#index sensi analysis



getsensindex <- function(param, type="arith")
{
	sapply(1:NROW(param), function(i)
		{
			par <- param[i, ]
			#print(par)
			xstar <- getNE(par)
			ifelse(type == "arith", mean(xstar), prod(sqrt(abs(xstar))))
		}
	)
	
}

X1 <- t(replicate(n, genpar()))
X2 <- t(replicate(n, genpar()))
Y1bar <- getsensindex(X1)
Y1geo <- getsensindex(X1, type="geo")

scr3 <- src(X1, Y1bar)
scr3
pcc3 <- pcc(data.frame(X1), Y1bar)


r3 <- sobol(getsensindex, X1, X2, order=1)

indexarith <- cbind(scr3$SRC, pcc3$PCC, r3$S)
colnames(indexarith) <- paste("Arith. Index.", c("SCR", "PCC", "Sobol"), sep="-")

scr4 <- src(X1, Y1geo)
pcc4 <- pcc(data.frame(X1), Y1bar)
r4 <- sobol(getsensindex, X1, X2, order=1, type="geo")

indexgeo <- cbind(scr4$SRC, pcc4$PCC, r4$S)
colnames(indexgeo) <- paste("Arith. Geo.", c("SCR", "PCC", "Sobol"), sep="-")

#--------------------------------
#function sensi analysis

basepar <- c(d1=20, d2=20, lam1=4, lam2=4, rho1=1, rho2=1)
xbase <- getNE(basepar)


getsensiNIF <- function(param)
{
	sapply(1:NROW(param), function(i)
		{
			par <- param[i, ]
			#print(par)
			xstar <- getNE(par)
			NIF(xbase, xstar, par)
		}
	)
	
}

X1 <- t(replicate(n, genpar()))
X2 <- t(replicate(n, genpar()))
Y1 <- getsensiNIF(X1)

scr5 <- src(X1, Y1)
scr5
pcc5 <- pcc(data.frame(X1), Y1)

r5 <- sobol(getsensiNIF, X1, X2, order=1)


funcNIF <- cbind(scr5$SRC, pcc5$PCC, r5$S)
colnames(funcNIF) <- paste("Func. NIF", c("SCR", "PCC", "Sobol"), sep="-")



getsensiSumobj <- function(param)
{
	sapply(1:NROW(param), function(i)
		{
			par <- param[i, ]
			#print(par)
			xstar <- getNE(par)
			o1 <- obj(xstar, 1, par["d1"], par["lam1"], par["rho1"]) 
			o2 <- obj(xstar, 2, par["d2"], par["lam2"], par["rho2"])			
			o1+o2
		}
	)
	
}

X1 <- t(replicate(n, genpar()))
X2 <- t(replicate(n, genpar()))
Y1 <- getsensiSumobj(X1)

scr6 <- src(X1, Y1)
pcc6 <- pcc(data.frame(X1), Y1)

r6 <- sobol(getsensiNIF, X1, X2, order=1)

funcSumObj <- cbind(scr6$SRC, pcc6$PCC, r6$S)
colnames(funcSumObj) <- paste("Func. Sum. Obj", c("SCR", "PCC", "Sobol"), sep="-")


resfinal <- cbind(playcomp, playobj, indexarith, indexgeo, funcNIF, funcSumObj)

paramname <- rownames(resfinal)

ranklabel <- function(x)
{
	x <- abs(x)
	paramname[order(x, decreasing=TRUE)]
}

paramranks <- apply(resfinal, 2, ranklabel)
rownames(paramranks) <- paste("top", 1:length(paramname))

paramranks[, 1:6]
paramranks[, 1:6+6]


