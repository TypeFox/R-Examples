### R code from vignette source 'GNE-howto.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: load
###################################################
library(GNE)


###################################################
### code chunk number 2: argphi
###################################################
myarg <- list(C=c(2, 3), D=c(4,0))
dimx <- c(1, 1)
#Gr_x_j O_i(x)
grobj <- function(x, i, j, arg)
{
	dij <- 1*(i == j)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[i] - arg$C[i])*(x[other] - arg$D[i])^4*dij 
	res + 4*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^3*(1-dij) 
}

dimlam <- c(1, 1)
#g_i(x)
g <- function(x, i)
	ifelse(i == 1, sum(x[1:2]) - 1, 2*x[1]+x[2]-2)
#Gr_x_j g_i(x)
grg <- function(x, i, j)
	ifelse(i == 1, 1, 1 + 1*(i == j))


###################################################
### code chunk number 3: argJacF
###################################################
#Gr_x_k Gr_x_j O_i(x)
heobj <- function(x, i, j, k, arg)
{
	dij <- 1*(i == j)
	dik <- 1*(i == k)
	other <- ifelse(i == 1, 2, 1)
	res <- 2*(x[other] - arg$D[i])^4*dij*dik 
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*dij*(1-dik)
	res <- res + 8*(x[i] - arg$C[i])*(x[other] - arg$D[i])^3*(1-dij)*dik
	res + 12*(x[i] - arg$C[i])^2*(x[other] - arg$D[i])^2*(1-dij)*(1-dik)
}
#Gr_x_k Gr_x_j g_i(x)
heg <- function(x, i, j, k) 0


###################################################
### code chunk number 4: testGNE
###################################################
set.seed(1234)
z0 <- rexp(sum(dimx)+sum(dimlam))
GNE.nseq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	compl=phiFB, gcompla=GrAphiFB, gcomplb=GrBphiFB, method="Newton", 
	control=list(trace=0))


###################################################
### code chunk number 5: trueGNE
###################################################
#list of true GNEs
trueGNE <- rbind(c(2, -2, 0, 5*2^5),
	c(-2, 3, 8, 0),
	c(0, 1, 4*3^4, 0),
	c(1, 0, 2^9, 6))
colnames(trueGNE) <- c("x1", "x2", "lam1", "lam2")
rownames(trueGNE) <- 1:4
print(trueGNE)


###################################################
### code chunk number 6: singjac
###################################################
z0 <- c(0, 0, 1, 1)
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiMin, gcomplb=GrBphiMin)


###################################################
### code chunk number 7: singjac2
###################################################
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiFB, gcomplb=GrBphiFB)
jacSSR(z0, dimx, dimlam, heobj=heobj, myarg, constr=g, grconstr=grg, 
	heconstr=heg, gcompla=GrAphiKK, gcomplb=GrBphiKK, argcompl=3/2)


###################################################
### code chunk number 8: testGNEceq
###################################################
z0 <- 1+rexp(sum(dimx)+2*sum(dimlam))
GNE.ceq(z0, dimx, dimlam, grobj=grobj, myarg, heobj=heobj, myarg, 
	constr=g, grconstr=grg, heconstr=heg, 
	method="PR", control=list(trace=0))


###################################################
### code chunk number 9: testNI
###################################################
#O_i(x)
obj <- function(x, i, arg)
  (x[i] - arg$C[i])^2*(x[-i] - arg$D[i])^4
#g(x)
gtot <- function(x)
  sum(x[1:2]) - 1
#Gr_x_j g(x)
jacgtot <- function(x)
	cbind(1, 1)

z0 <- rexp(sum(dimx))

GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, 
         jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), 
         problem="NIR", merit="NI")


GNE.fpeq(z0, dimx, obj, myarg, grobj, myarg, heobj, myarg, gtot, NULL, 
         jacgtot, NULL, silent=TRUE, control.outer=list(maxit=10), 
         problem="VIR", merit="VI")


