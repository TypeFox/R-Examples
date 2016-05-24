## test column order options of oa.design
require(DoE.base)

P3.3(oa.design(L18, nlevels=c(2,3,3,3),columns="order"), detailed=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3),columns="min3"), detailed=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3),columns="min34"), detailed=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3),columns="minRPFT"), detailed=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3),columns="minRelProjAberr"), detailed=TRUE)


P3.3(oa.design(L18, nlevels=c(2,3,3,3,3),columns="min3"), parft=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3),columns="min34"), parft=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3),columns="minRPFT"), parft=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3),columns="minRelProjAberr"), parft=TRUE)

P3.3(oa.design(L18, nlevels=c(3,3,2,3,3),columns="min3"), parftdf=TRUE)
P3.3(oa.design(L18, nlevels=c(3,3,2,3,3),columns="min34"), parftdf=TRUE)
P3.3(oa.design(L18, nlevels=c(3,3,2,3,3),columns="minRPFT"), parftdf=TRUE)
P3.3(oa.design(L18, nlevels=c(3,3,2,3,3),columns="minRelProjAberr"), parftdf=TRUE)

P3.3(oa.design(L18, nlevels=c(2,3,3,3,3,3),columns="min3"), rela=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3,3),columns="min34"), rela=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3,3),columns="minRPFT"), rela=TRUE)
P3.3(oa.design(L18, nlevels=c(2,3,3,3,3,3),columns="minRelProjAberr"), rela=TRUE)

GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)))
GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)), arft=FALSE)
GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)), scft=FALSE)
GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)), arft=FALSE, scft=FALSE)
GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)), arft=FALSE, scft=FALSE, cancor=TRUE)
GRind(oa.design(L18, nlevels=c(2,3,3,3,3,3,3)), cancor=TRUE)

## interesting, but take too long
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="min3"))
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="min3.rela"), rela=TRUE)
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="min34"))
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="min34.rela"), rela=TRUE)
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="minRPFT"), rela=TRUE)
#P3.3(oa.design(L36.2.11.3.12, nlevels=c(2,2,2,3,3,3),columns="minRelProjAberr"), rela=TRUE)

## also interesting but also take too long
#P3.3(oa.design(L32.2.10.4.7, nlevels=c(2,2,2,4,4,4,4,4), columns="order"), rela=TRUE)
#P3.3(oa.design(L32.2.10.4.7, nlevels=c(2,2,2,4,4,4,4,4), columns="min34"), rela=TRUE)
#P3.3(oa.design(L32.2.10.4.7, nlevels=c(2,2,2,4,4,4,4,4), columns="min34.rela"), rela=TRUE)
#P3.3(oa.design(L32.2.10.4.7, nlevels=c(2,2,2,4,4,4,4,4), columns="minRPFT"), rela=TRUE)
#P3.3(oa.design(L32.2.10.4.7, nlevels=c(2,2,2,4,4,4,4,4), columns="minRelProjAberr"), rela=TRUE)