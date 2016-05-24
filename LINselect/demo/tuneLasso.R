library("LINselect")

# simulate data with
ex <- simulData(p=100,n=100,r=0.8,rSN=5)

ex1.tuneLasso <- tuneLasso(ex$Y,ex$X)

data(diabetes)
attach(diabetes)
ex.diab <- tuneLasso(y,x2)
detach(diabetes)



