## ----install, eval=FALSE-------------------------------------------------
#  install.packages('nat.nblast', dependencies = TRUE)

## ------------------------------------------------------------------------
library(nat.nblast)
rgl::setupKnitr()

## ---- message=FALSE------------------------------------------------------
library(nat)
kcscores <- nblast_allbyall(kcs20)

## ------------------------------------------------------------------------
hckcs <- nhclust(scoremat=kcscores)
library(dendroextras)
dkcs <- colour_clusters(hckcs, k=3)

## ------------------------------------------------------------------------
labels(dkcs) <- with(kcs20[labels(dkcs)], type)
plot(dkcs)

## ---- rgl=TRUE, message=FALSE--------------------------------------------
plot3d(hckcs, k=3, db=kcs20, soma=T)
par3d(userMatrix=diag(c(1,-1,-1,1), 4))
plot3d(MBL.surf, alpha=.1)

