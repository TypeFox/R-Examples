if(!require("GNE"))stop("this test requires package GNE.")

#-------------------------------------------------------------------------------
# (4) Example of GNE with 4 solutions(!)
#-------------------------------------------------------------------------------

myarg <- list(C=c(2, 3), D=c(4,0))

dimx <- c(1, 1)


#O_i(x)
obj <- function(x, i, arg)
	(x[i] - arg$C[i])^2*(x[-i] - arg$D[i])^4




#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	dij <- 1*(i == j)
	other <- ifelse(i == 1, 2, 1)
	2*(x[i] - arg$C[i])*(x[other] - arg$D[i])^4*dij + 4*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^3*(1-dij) 
}
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
	dij <- 1*(i == j)
	dik <- 1*(i == k)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[other] - arg$D[i])^4*dij*dik + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*dij*(1-dik)
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*(1-dij)*dik
	res + 12*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^2*(1-dij)*(1-dik)
}

dimlam <- c(1, 1)
#g(x)
gtot <- function(x)
	sum(x[1:2]) - 1
#	c(sum(x[1:2]) - 1, 2*x[1]+x[2]-2)
#Gr_x_j g(x)
jacgtot <- function(x)
	cbind(1, 1)
#	cbind(c(1, 1), c(2, 1))




z0 <- rexp(sum(dimx))

fpNIR(z0, dimx, obj, myarg, gtot, NULL, grobj, myarg, jacgtot, NULL)

GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), problem="NIR", merit="NI")

z0 <- rexp(sum(dimx))

fpVIR(z0, dimx, obj, myarg, gtot, NULL, grobj, myarg, jacgtot, NULL)

GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), problem="VIR", merit="VI")

