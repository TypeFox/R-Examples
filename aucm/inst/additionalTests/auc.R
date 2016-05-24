
seed=26
seed=16
seed=3
dat.train = sim.dat.1(n=200, seed=seed, add.outliers=TRUE) 
fits=list()
fits[[1]]=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="Lin")
fits[[2]]=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="MH")
fits[[3]]=sauc.phi(y~x1+x2, dat.train,constrain.method="beta1",h.method="Lin")
fits[[4]]=sauc.phi(y~x1+x2, dat.train,constrain.method="beta1",h.method="MH") # not a good combination of constrain.method and h.method
sapply(fits, function(x) ratio(x)[2])

# compared to rlogit
fit.rlogit=BYlogreg(y~x1+x2, dat.train)
coef(fit.rlogit)[3]/coef(fit.rlogit)[2]

# compared to rauc
fit1=rauc.linear (y~x1+x2, dat.train, s=1,lambda=1,beta.init=c(1,1),minQuad.tol=1e-3)
ratio(fit1)
fit2=rauc.linear (y~x1+x2, dat.train, s=1,lambda=1,beta.init=coef(fit.rlogit)[-1],minQuad.tol=1e-3)
ratio(fit2)

# compare Shuxin and Ying's code
dat.train = sim.dat.1(n=200, seed=26, add.outliers=TRUE) 
fit1=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="Vexler",start.method="rlogit")
#fit2=calLin.sy(dat.train,start.method="rlogistic") # calLin.sy not found anymore

# explosion
seed=954
dat.train = sim.dat.1(n=200, seed=seed, add.outliers=TRUE) 
fit.1=sauc.phi(y~x1+x2, dat.train,constrain.method="L2",h.method="Lin")
ratio(fit.1)
