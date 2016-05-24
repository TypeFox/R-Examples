### R code from vignette source 'Elston-Stewart.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Elston-Stewart.Rnw:22-23
###################################################
options(continue=" ", prompt = " ", SweaveHooks=list(fig=function() par(mar=c(5.1,4.1,3.1,2.1))), width=90)


###################################################
### code chunk number 2: prompton
###################################################
options(prompt="> ", continue = "+ ");


###################################################
### code chunk number 3: promptoff
###################################################
options(prompt=" ", continue=" ");


###################################################
### code chunk number 4: Elston-Stewart.Rnw:34-35
###################################################
options(prompt="> ", continue = "+ ");


###################################################
### code chunk number 5: desc
###################################################
require(ElstonStewart)
options(width = 90)
desc <- packageDescription("ElstonStewart")


###################################################
### code chunk number 6: Elston-Stewart.Rnw:94-95
###################################################
modele.di$proba.g


###################################################
### code chunk number 7: Elston-Stewart.Rnw:99-100
###################################################
modele.di$trans


###################################################
### code chunk number 8: Elston-Stewart.Rnw:105-106
###################################################
modele.di$p.pheno


###################################################
### code chunk number 9: Elston-Stewart.Rnw:122-124
###################################################
data(conrad2)
conrad2


###################################################
### code chunk number 10: Elston-Stewart.Rnw:130-135
###################################################
genotypes <- c( rep(list(0:2), 21), 2 )

X <- es.pedigree( id = conrad2$id, father = conrad2$father, mother = conrad2$mother, 
      sex = conrad2$sex, pheno = rep(0, 22), geno = genotypes )
X


###################################################
### code chunk number 11: Elston-Stewart.Rnw:140-141
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(X)


###################################################
### code chunk number 12: Elston-Stewart.Rnw:146-148
###################################################
r <- Elston(X, modele.di, list(p = 0.98))
r$result


###################################################
### code chunk number 13: Elston-Stewart.Rnw:154-158
###################################################
# using the memoization...
system.time(r <- Elston(X, modele.di, list(p = 0.98)))
system.time(r <- Elston(X, modele.di, list(p = 0.98), r$mem))
system.time(r <- Elston(X, modele.di, list(p = 0.99), r$mem))


###################################################
### code chunk number 14: Elston-Stewart.Rnw:164-169
###################################################
modele.rec <- list( name = "recessive", proba.g = modele.di$proba.g, 
   trans = modele.di$trans,
   p.pheno = function(x, g, theta) 
       ifelse( is.na(x) | (x == 1 & g == 2) | (x == 0 & g < 2) , 1, 0) 
   )


###################################################
### code chunk number 15: Elston-Stewart.Rnw:173-179
###################################################
genotypes <- rep(list(0:2), 22)
X <- es.pedigree( id = conrad2$id, father = conrad2$father, mother = conrad2$mother, 
      sex = conrad2$sex, pheno = c( rep(NA, 21), 1), geno = genotypes )

r <- Elston(X, modele.rec, list(p = 0.98), r$mem)
r$result


###################################################
### code chunk number 16: Elston-Stewart.Rnw:184-189
###################################################
X <- es.pedigree( id = conrad2$id, father = conrad2$father, mother = conrad2$mother,
      sex = conrad2$sex, pheno = c( rep(0, 21), 1), geno = genotypes )

r <- Elston(X, modele.rec, list(p = 0.98), r$mem)
r$result


###################################################
### code chunk number 17: Elston-Stewart.Rnw:194-199
###################################################
genotypes <- c( rep(list(0:1), 21), 2 )
X <- es.pedigree( id = conrad2$id, father = conrad2$father, mother = conrad2$mother,
      sex = conrad2$sex, pheno = rep(0, 22), geno = genotypes )
r <- Elston(X, modele.di, list(p = 0.98), r$mem)
r$result


###################################################
### code chunk number 18: Elston-Stewart.Rnw:206-208
###################################################
data(fams)
head(fams,15)


###################################################
### code chunk number 19: Elston-Stewart.Rnw:213-227
###################################################
fam.ids <- unique(fams$fam);

# creating a list of genotypes corresponding to individuals in fam.ids
# genotype is NA -> 0, 1 or 2
genotypes <- lapply( fams$genotype, function(x) if(is.na(x)) 0:2 else x )

X <- vector("list", length(fam.ids))
for(i in seq_along(fam.ids))
{
  w <- which(fams$fam == fam.ids[i])
  X[[i]] <- es.pedigree( id = fams$id[w], father = fams$father[w],
      mother = fams$mother[w], sex = fams$sex[w], pheno = rep(0, length(w)), 
      geno = genotypes[w], famid = fam.ids[i] )
}


###################################################
### code chunk number 20: Elston-Stewart.Rnw:231-233
###################################################
# computing the log-likelihood for a single value p
Likelihood(X, modele.di, theta = list( p=0.5), n.cores=1 )


###################################################
### code chunk number 21: Elston-Stewart.Rnw:239-243
###################################################
getOption("SweaveHooks")[["fig"]]()
# computing the likelihood for a vector p
p <- seq(0,1,length=501)
L <- Likelihood(X, modele.di, theta = list( p=p ), n.cores=1 ) 
plot( p, exp(L), type="l")


###################################################
### code chunk number 22: Elston-Stewart.Rnw:249-253
###################################################
# running an optimization algorithm
# Elston-Stewart is ran several times
# here we run the algorithm with 2 cores
optimize( function(p) -Likelihood(X, modele.di, theta = list( p=p ), n.cores=2 ) , c(0.35,0.45) )


###################################################
### code chunk number 23: Elston-Stewart.Rnw:258-259
###################################################
es.stopCluster()


