### R code from vignette source 'pedgene.Rnw'

###################################################
### code chunk number 1: desc
###################################################
options(width = 100)
desc <- packageDescription("pedgene")


###################################################
### code chunk number 2: pedgene.Rnw:43-44 (eval = FALSE)
###################################################
## library(pedgene)


###################################################
### code chunk number 3: pedgene.Rnw:47-49
###################################################
.libPaths(new = c("/projects/bsi/gentools/R/lib301", .libPaths()))
library(pedgene)


###################################################
### code chunk number 4: exdata
###################################################

data(example.ped)
example.ped[c(1:10,31:39),]
data(example.geno)
dim(example.geno)
example.geno[,c(1:2,10:14)]
data(example.map)
example.map


###################################################
### code chunk number 5: ccbase
###################################################
pg.cc <- pedgene(ped=example.ped, geno=example.geno, map=example.map)

print(pg.cc)


###################################################
### code chunk number 6: cccov
###################################################
binfit <-  glm(trait ~ (sex-1),data=example.ped, family="binomial")

example.ped$trait.adjusted[!is.na(example.ped$trait)] <- fitted.values(binfit) 
example.ped[1:10,]

pg.cc.adj <- pedgene(ped=example.ped, geno=example.geno, map=example.map)

summary(pg.cc.adj, digits=4)


###################################################
### code chunk number 7: continuous
###################################################
set.seed(10)
cont.ped <- example.ped[,c("ped", "person", "father", "mother", "sex")]
beta.sex <- 0.3
cont.ped$trait <- (cont.ped$sex-1)*beta.sex + rnorm(nrow(cont.ped))

pg.cont <- pedgene(ped = cont.ped, geno = example.geno, map = example.map)
print(pg.cont, digits=4)


###################################################
### code chunk number 8: continuousadj
###################################################
gausfit <-  glm(trait ~ (sex-1),data=cont.ped, family="gaussian")

cont.ped$trait.adjusted <- fitted.values(gausfit) 
cont.ped[1:10,]

pg.cont.adj <- pedgene(ped=cont.ped, geno=example.geno, map=example.map)

summary(pg.cont.adj, digits=5)


###################################################
### code chunk number 9: pedgene.Rnw:184-185
###################################################
toLatex(sessionInfo())


