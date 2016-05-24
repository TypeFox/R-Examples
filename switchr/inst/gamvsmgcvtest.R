##http://stackoverflow.com/questions/22126611/r-package-conflict-between-gam-and-mgcv

library(switchr)
                                        # fresh session
t.s1 <- search()
t.lN1 <- loadedNamespaces()

# some dummy data
data <-data.frame(is.exc=sample(x=c(0,1),size=100,replace=T),
year=1:100,doy=rep(1:5,times=20))
t.dof <- 2

# everything works fine
library(gam)
t.gam1 <- gam::gam(is.exc~s(year,df=t.dof)+s(doy,df=t.dof),data=data,family=poisson)
t.pred1 <- gam::predict.gam(t.gam1,newdata=data,type='terms')
detach('package:gam',unload=T,character.only=T)
detach('package:splines',unload=T,character.only=T)

# compare attached packages and namespaces with fresh session
t.s2 <- search()
t.lN2 <- loadedNamespaces()
identical(t.s1,t.s2)
identical(t.lN1,t.lN2)
flushSession()

# attach and detach mgcv
library(mgcv)
detach('package:mgcv',unload=T,character.only=T)
unloadNamespace('nlme')
unloadNamespace('Matrix')
unloadNamespace('lattice')
unloadNamespace('grid')

# compare again attached packages and namespaces with fresh session
t.s2 <- search()
t.lN2 <- loadedNamespaces()
identical(t.s1,t.s2)
identical(t.lN1,t.lN2)

flushSession()

# use package gam again and produce errors
library(gam)
t.gam2 <- gam::gam(is.exc~s(year,df=t.dof)+s(doy,df=t.dof),data=data,family=poisson)
gam::summary.gam(t.gam2)
t.pred2 <- gam::predict.gam(t.gam2,newdata=data,type='terms')

# why do we have mgcv and friends in the namespace?
loadedNamespaces()
