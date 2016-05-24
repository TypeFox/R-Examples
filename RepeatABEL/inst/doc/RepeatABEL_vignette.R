### R code from vignette source 'RepeatABEL_vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: RepeatABEL_vignette.Rnw:106-108
###################################################
options(width=60)
options(continue="    ")


###################################################
### code chunk number 2: RepeatABEL_vignette.Rnw:148-153
###################################################
library(RepeatABEL)
#GenABEL object including IDs and marker genotypes
data(gen.data) 
#Phenotype data with repeated observations
data(Phen.Data) 


###################################################
### code chunk number 3: RepeatABEL_vignette.Rnw:158-160
###################################################
GWAS1 <- rGLS(y ~ age + sex, genabel.data = gen.data, 
phenotype.data = Phen.Data) 


###################################################
### code chunk number 4: RepeatABEL_vignette.Rnw:168-169
###################################################
summary(GWAS1)


###################################################
### code chunk number 5: RepeatABEL_vignette.Rnw:176-180
###################################################
est.hglm <- GWAS1@call$hglm
cat("Genotypic and permanent env. variance components:","\n", 
    est.hglm$varRanef,", resp.","\n",
    "The residual variance is", est.hglm$varFix,".","\n")


###################################################
### code chunk number 6: RepeatABEL_vignette.Rnw:186-191
###################################################
logVCE <- est.hglm$SummVC2[[1]]
cat("Estimate and SE for the genotypic variance on a natural log scale:", 
    "\n", logVCE[1],"and",logVCE[2],", resp.", "\n \n",
    "Confidence interval: [", exp( logVCE[1] - 1.96*logVCE[2] ) , "," , 
    exp( logVCE[1] + 1.96*logVCE[2] ) , "]" , "\n")


###################################################
### code chunk number 7: RepeatABEL_vignette.Rnw:197-200
###################################################
cat("Heritability:",
    est.hglm$varRanef[1]/(est.hglm$varFix + est.hglm$varRanef[2]), 
    "\n")


###################################################
### code chunk number 8: RepeatABEL_vignette.Rnw:211-220
###################################################
#The same results can be computed using the preFitModel as follows
fixed=y ~ age + sex
Mod1 <- preFitModel(fixed, random=~1|id, 
  genabel.data = gen.data, 
  phenotype.data = Phen.Data, 
  corStruc=list( id=list("GRM","Ind") )) 
GWAS1b <- rGLS(fixed, genabel.data = gen.data, 
phenotype.data = Phen.Data, V = Mod1$V) 
summary(GWAS1b)


###################################################
### code chunk number 9: RepeatABEL_vignette.Rnw:229-239 (eval = FALSE)
###################################################
## # In this example there are 6 nests and 60 observations per nest
## Phen.Data$nest <- rep(1:6, each=60)
## 
## # A model including polygenic effects,  
## # permanent environmental effects, and nest effect as random
## Mod2 <- preFitModel(fixed, random=~1|id + 1|nest, 
## genabel.data = gen.data, phenotype.data = Phen.Data, 
## corStruc=list( id=list("GRM","Ind") , nest=list("Ind")) )
## GWAS2 <- rGLS(fixed, genabel.data = gen.data, 
## phenotype.data = Phen.Data, V = Mod2$V) 


###################################################
### code chunk number 10: RepeatABEL_vignette.Rnw:246-253
###################################################
D= matrix(0,6,6)
D[1,2] = D[2,1] = 1
D[5,6] = D[6,5] = 1
D[2,4] = D[4,2] = 1
D[3,5] = D[5,3] = 1
D[1,6] = D[6,1] = 1
D[3,4] = D[4,3] = 1


###################################################
### code chunk number 11: RepeatABEL_vignette.Rnw:258-264 (eval = FALSE)
###################################################
## Mod3 <- preFitModel(y ~ age + sex, random=~1|id + 1|nest, 
## genabel.data = gen.data, phenotype.data = Phen.Data, 
## corStruc=list( id=list("GRM","Ind") , 
## nest=list("CAR")), Neighbor.Matrix=D )
## GWAS2b <- rGLS(fixed, genabel.data = gen.data, 
## phenotype.data = Phen.Data, V = Mod3$V) 


###################################################
### code chunk number 12: RepeatABEL_vignette.Rnw:272-274
###################################################
VC.poly <- VC.perm <- VC.res <- 1
n.obs <- rep(4, nids(gen.data)) 


###################################################
### code chunk number 13: RepeatABEL_vignette.Rnw:279-283 (eval = FALSE)
###################################################
## Phen.Sim <- simulate_PhenData(y ~ 1, 
##   genabel.data = gen.data, 
##   n.obs = n.obs, SNP.eff = 2, SNP.nr = 1000, 
##   VC = c(VC.poly, VC.perm,VC.res))


###################################################
### code chunk number 14: RepeatABEL_vignette.Rnw:288-290 (eval = FALSE)
###################################################
## GWAS.sim1 <- rGLS(y ~ 1, genabel.data = gen.data, 
## phenotype.data = Phen.Sim)


###################################################
### code chunk number 15: RepeatABEL_vignette.Rnw:295-298 (eval = FALSE)
###################################################
## Phen.Sim <- simulate_PhenData(y ~ sex, 
## genabel.data = gen.data, n.obs = n.obs, SNP.eff = 2, 
## SNP.nr = 1000, VC = c(VC.poly, VC.perm, VC.res), beta = c(0,1))


###################################################
### code chunk number 16: RepeatABEL_vignette.Rnw:303-306 (eval = FALSE)
###################################################
## GWAS.sim1 <- rGLS(y ~ sex, genabel.data = gen.data, 
## phenotype.data = Phen.Sim)
## plot(GWAS.sim1, main="Simulation results")


###################################################
### code chunk number 17: RepeatABEL_vignette.Rnw:317-345 (eval = FALSE)
###################################################
## Phen.Sim <- simulate_PhenData(y ~ sex, 
## genabel.data = gen.data, n.obs = n.obs, SNP.eff = 2, 
## SNP.nr = 1000, VC = c(VC.poly, VC.perm, VC.res), beta = c(0,1))
## ####### PRODUCE YEAR EFFECTS FOR EACH INDIVIDUAL
## sd.year = 1 #Standard Deviation of Year Effects
## beta <- rnorm(max(n.obs), 0, sd.year) #Simulated Year Effects
## year.effects <- years <- NULL
## for (i in 1:length(n.obs)) {
##   yr.i <- sort(sample(1:max(n.obs), n.obs[i]))
##   years <- c(years, yr.i)
## 	year.effects <- c(year.effects, beta[yr.i])
## }
## ########################
## #### A FUNCTION TO ADD A VARIABLE TO A LIST
## add.var <- function(x, add.new, new.name) {
##   x[[length(x) + 1]] <- add.new
## 	names(x)[length(x)] <- new.name
## 	return(x)
## }
## #######################
## #ADDS A NEW PHENOTYPE WITH YEAR EFFECTS ADDED
## Phen.Sim <- add.var(Phen.Sim, Phen.Sim$y + year.effects, "y.yrs")  
## #ADD YEAR AS FACTOR
## Phen.Sim <- add.var(Phen.Sim, as.factor(years), "Years") 
## #######################
## #RUN THE ANALYSIS
## GWAS.sim1 <- rGLS(y.yrs ~ sex + Years, genabel.data = gen.data, 
## phenotype.data = Phen.Sim)


###################################################
### code chunk number 18: RepeatABEL_vignette.Rnw:352-359 (eval = FALSE)
###################################################
## Phen.Sim <- simulate_PhenData(y ~ sex, 
## genabel.data = gen.data, n.obs = n.obs, SNP.eff = 2, 
## SNP.nr = 1000, VC = c(VC.poly, VC.perm, VC.res), beta = c(0,1))
## tau = 0
## Phen.Sim$y_observed <- as.numeric(Phen.Sim$y > tau)
## GWAS.sim1 <- rGLS(y_observed ~ sex, genabel.data = gen.data, 
## phenotype.data = Phen.Sim)


