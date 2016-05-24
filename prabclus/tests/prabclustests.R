library(prabclus)
options(digits=4)

data(kykladspecreg)
data(nb)
set.seed(1234)
x <- prabinit(prabmatrix=kykladspecreg, neighborhood=nb)

p1 <- prabtest(x, times=3, pd=0.35, ignore.richness=TRUE)
p2 <- prabtest(x, times=3, pd=0.35, teststat="lcomponent")
p3 <- prabtest(x, times=3, pd=0.35, teststat="isovertice")
p4 <- prabtest(x, times=3, pd=0.35, teststat="nn", sf.sim=TRUE)
p5 <- prabtest(x, times=3, pd=0.35, teststat="inclusions")
summary(p1)
summary(p2)
summary(p3)
summary(p4)
summary(p5)

data(veronica)
vnb <- coord2dist(coordmatrix=veronica.coord[1:50,], cut=20,
                  file.format="decimal2",neighbors=TRUE)
vei <- prabinit(prabmatrix=veronica[1:50,],
                neighborhood=vnb$nblist,nbbetweenregions=FALSE,
                distance="jaccard")
print(vei)

library(spdep)
data(siskiyou)
x <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb,
            distance="logkulczynski")
build.nblist(x)
a1 <- abundtest(x, times=5, p.nb=0.0465)
a2 <- abundtest(x, times=5, p.nb=0.0465, teststat="groups",
                groupvector=siskiyou.groups)
# These settings are chosen to make the example execution
# faster; usually you will use abundtest(x).
summary(a1)
summary(a2)

options(digits=2)
prab.sarestimate(x)
regpop.sar(x, p.nb=0.046)
options(digits=4)

x <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb, distance="none",toprab=TRUE,toprabp=0.5)
x2 <- prabinit(prabmatrix=siskiyou, neighborhood=siskiyou.nb, distance="none",toprab=TRUE,toprabp=0)
x$prab
x2$prab


