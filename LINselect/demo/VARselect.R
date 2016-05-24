library("LINselect")

# simulate data with
# beta=c(rep(2.5,5),rep(1.5,5),rep(0.5,5),rep(0,p-15))
ex <- simulData(p=100,n=100,r=0.8,rSN=5)

ex1.Vselect <- VARselect(ex$Y,ex$X,exhaustive.dmax=2)

data(diabetes)
attach(diabetes)
ex.diab <- VARselect(y,x2,exhaustive.dmax=5)
detach(diabetes)



