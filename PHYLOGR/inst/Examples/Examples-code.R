#### Examples-code.R   original: 30-June-2000;
#### a few later (August, October, November, March-2001) modifications.



### This file contains all the commands used to generate the examples
### shown in Díaz-Uriarte & Garland, in review.

### All the files that are read are in the 'Examples' directory
### of the package (which is also where this file you are reading
### was originally located).

### For reasons of space, the simulated data files (*.sim) are small
### files with only 100 simulated data sets (in contrast to the 5,000
### simulations used in the paper).  To obtain more realistic
### results, you might want to run PDSIMUL and generate a larger
### number of simulations. (For completeness: the seeds we used
### were 1000 for 49lbr.inp, 100 for 49ms.pdi and 20 for 49hmt.pdi).


### To reproduce the results of the paper, the simplest would be
### to either copy all these files to a new directory and launch
### R from there, or launch R from the 'Examples' directory (this
### might not work under Linux if you don't have writing permissions
### to that directory and want to generate postscript files without
### specifying another directory).

### You can run this file either one command at a time (cutting and
### pasting or, better yet, using ESS), or you can directly source the
### file into R.
### To use "source" do:
### > source("Examples-code.R",echo=TRUE,print.eval=TRUE,max.deparse.length=500)
### that will run the file, echo the commands, and return output.



options(digits=4) # number of output digits
library(PHYLOGR) # load our package



###########################
###########################



ANCOVA.sim <- read.sim.data("49lbr.sim", "49lbr.inp", variable.names = 	c("body.mass","home.range"), other.variables = data.frame(clade=c(rep("Carnivore",19),rep("Herbivore",30)))) # this reads in the simulated data file, and adds a column with the herbivore vs. carnivore indicator

ancova.fit <- lm.phylog(home.range ~ body.mass + clade, data = ANCOVA.sim) # fit the model

summary(ancova.fit)

#############################
##############################

garland.janis.ic <- cbind(read.table("49ms.fic")[,c(3,4)], read.table("49hmt.fic")[,c(3,4)]) # read fic files

branch.lengths <- read.table("49ms.fic")[,5] # the branch lengths

garland.janis.ic <- garland.janis.ic/branch.lengths # standardize contrasts

names(garland.janis.ic) <- c("body.mass", "running.speed", "hind.l.length", "mtf.ratio") # name the variables

garland.janis.ic[garland.janis.ic$body.mass<0,] <- -1 * garland.janis.ic[ garland.janis.ic$body.mass<0, ] # to positivize contrasts

garland.janis.ic$branch.lengths <- branch.lengths # add branch lengths 	to the file, for completennes

garland.janis.ic$names.contr <- as.factor(read.table("49ms.fic")[,1])

garland.janis.ic$clade.contr <- as.factor( c("root",rep("Carnivore",18), rep("Herbivore",29))) # create indicator for clade


attach(garland.janis.ic)

# correlation of residuals from regressions on body mass
resid.hll <- resid(lm(hind.l.length ~ body.mass - 1)) # this model formula uses a regression through the origin (the -1).

resid.speed <- resid(lm(running.speed ~ body.mass - 1))


rho1 <- cor.origin(resid.hll,resid.speed)
rho1 # the correlation through the origin
2*(1-pt(sqrt(46)*rho1/sqrt(1-rho1^2),46)) # the p-value, using the standard formula but providing the correct df's.


# the same, but using multiple regression
plot(body.mass,running.speed)
plot(body.mass,hind.l.length)
plot(hind.l.length,running.speed)
fit1 <- lm(running.speed ~ body.mass + hind.l.length - 1) # note also 	regression through the origin
plot(fit1)
summary(fit1)


# Excluding the polar bear-grizzly contrast
# using residuals:
resid.hll.nbear <- resid(lm(hind.l.length ~ body.mass - 1, subset = 	names.contr!="Tm-Ur"))

resid.speed.nbear <- resid(lm(running.speed ~ body.mass - 1, subset=names.contr!="Tm-Ur"))

rho2 <- cor.origin(resid.hll.nbear,resid.speed.nbear)
rho2 # the correlation through the origin
2*(1-pt(sqrt(45)*rho2/sqrt(1-rho2^2),45)) # the p-value, using the standard formula but providing the correct df's.



# with multiple regression
fit2 <- lm(running.speed ~ body.mass + hind.l.length - 1, subset=names.contr!="Tm-Ur")
plot(fit2)
summary(fit2)

# is MT/F relevant after including hindl.lenght? remember: need to look only at coefficient for mtf.
summary(lm(running.speed ~ body.mass + hind.l.length + mtf.ratio - 1, subset=names.contr!="Tm-Ur"))


anova(lm(running.speed ~ body.mass + hind.l.length + mtf.ratio - 1, subset=names.contr!="Tm-Ur")) # anova summary




GarlandJanis.sim <- read.sim.data(c("49ms.sim","49hmt.sim"), pdi.files = c("49ms.pdi","49hmt.pdi"), variable.names = c("body.mass","running.speed","hind.l.length","mtf.ratio"), other.variables = data.frame(clade = c(rep("Carnivore",19),rep("Herbivore",30)))) 

plot(lm(running.speed ~ body.mass + hind.l.length, data = GarlandJanis.sim, subset= sim.counter==0)) # regression residual plots

speed.hll <- lm.phylog(running.speed ~ body.mass + hind.l.length, data = GarlandJanis.sim)

summary(speed.hll) 


#postscript(file="Figure1.ps") # to generate a postscript file; you won't see the figure in the screen.
par(mfrow=c(2,3)) # 2x3 figures in the page, so that all fit in one 	indow
par(cex.lab=1.5) # to make labels larger
plot(speed.hll) # plot the analyses from lm.phylog
#dev.off()

speed.hll.nbear <- lm.phylog(running.speed ~ body.mass + hind.l.length, data=GarlandJanis.sim, exclude.tips=c("Tm","Ur"))

summary(speed.hll.nbear)



speed.hll.nbear.int <- lm.phylog(running.speed ~ body.mass + hind.l.length*clade, data=GarlandJanis.sim, exclude.tips=c("Tm","Ur"))

summary(speed.hll.nbear.int)



attach(GarlandJanis.sim) # to allow subseting of Tips for excluding 	either Carnivores or Herbivores

speed.mass.length.carn <- lm.phylog(running.speed ~ body.mass + hind.l.length, data=GarlandJanis.sim, exclude.tips= c(as.character(Tips[clade=="Herbivore" & sim.counter==0]), "Tm","Ur")) # only carnivores, except polar bear and grizzly

summary(speed.mass.length.carn)



# Testing whether MT/F ratio adds anything
speed.mtf <- lm.phylog(running.speed ~ body.mass + hind.l.length + mtf.ratio, data=GarlandJanis.sim, exclude.tips=c("Tm","Ur"))
summary(speed.mtf) # we are only interested here in the F value of mtf.ratio


############### PCA ###############

library(PHYLOGR) # load the package, if you haven't already done so.

LacertidSim <- read.sim.data(c("ifsmi.sim", "ihshw.sim", "iclag.sim", "icfxx.sim"), pdi.files=c("ifsmi.pdi", "ihshw.pdi", "iclag.pdi", "icfxx.pdi"), variable.names= c("svl","svl.matur", "hatsvl", "hatweight", "clutch.size", "age.mat","cl.freq", "xx"))

LacertidSim <- LacertidSim[,-10] # exclude last column

LacertidPCA <- prcomp.phylog(LacertidSim) 

summary(LacertidPCA)

#par(mfrow=c(2,4))
#par(cex.lab=1.2)
#plot(LacertidPCA) # shown in Figure 2

## or we can use:

#postscript(file="Figure2.ps") 
par(mfrow=c(2,4)) 
par(cex.lab=1.2) 
plot(LacertidPCA)
#dev.off() 




LacertidIC <- cbind(read.table("ifsmi.fic")[,c(3,4)], read.table("ihshw.fic")[,c(3,4)], read.table("iclag.fic")[,c(3,4)], read.table("icfxx.fic")[,3]) # the (unstandardized) contrasts

stand <- read.table("ifsmi.fic")[,5] # the square root of the sum of 	branch lengths

LacertidIC <- LacertidIC/stand # compute the standardized contrasts

LacertidIC$contr <- read.table("ifsmi.fic")[,1] # name of the contrast

names(LacertidIC) <- c("svl","svl.matur", "hatsvl", "hatweight", "clutch.size", "age.mat","cl.freq","ICcontr")

cor.for.ic.pca <- matrix(nrow=7,ncol=7) 

for (i in 1:7) for (j in 1:7) cor.for.ic.pca[i,j] <- cor.origin(LacertidIC[[i]],LacertidIC[[j]]) # to obtain the 	correlation matrix based on regressions through the origin

ic.pca <- svd(cor.for.ic.pca) #so that the eigenvectors are normalized

cor.with.factors <- t(sqrt(ic.pca$d) * t(ic.pca$u))

cor.with.factors <- as.data.frame(cor.with.factors)

ic.pca$d # these are the eigenvalues


100*ic.pca$d/length(ic.pca$d) # percentage of variance


names(cor.with.factors) <- paste("PC",seq(1:7),sep="") # nicer names for output

row.names(cor.with.factors)<-names(LacertidIC)[-8] 

cor.with.factors # the correlation between the variables and the components


######################################
#############     ####################
############# GLS ####################
#############     ####################
######################################


#### Lacertids

LacertidOriginal <- LacertidSim[LacertidSim$sim.counter==0,-1] # we 	could read the pdi data again, but this is another simple way, 	since the original data are in the file with the simulated data

Lacertid.varcov <- read.phylog.matrix("ifsmi.dsc") # from pddist, option 5, in matrix form, with header

# attach the data set, so variables can be accessed directly 
# without reference to the data set

detach(GarlandJanis.sim)
attach(LacertidOriginal) 

t(apply(LacertidOriginal[,-c(1,2)],2,function(y){summary(phylog.gls.fit(svl,y,Lacertid.varcov))[[4]][2,c(1,2)]}))

#### Garland & Janis

Garland.Janis.Orig <- read.pdi.data(c("49ms.pdi","49hmt.pdi"), variable.names = c("body.mass", "running.speed", "hind.l.length", "mtf.ratio"))
Garland.Janis.Orig$clade <- as.factor(c(rep("Carnivore",19), rep("Herbivore",30)))
Garland.Janis.Cov <- read.phylog.matrix("49ms.dsc")

detach(LacertidOriginal)
detach(garland.janis.ic)
attach(Garland.Janis.Orig)

# all data
fit.gls.gj <- phylog.gls.fit(cbind(body.mass,hind.l.length), running.speed, Garland.Janis.Cov) # the model fitting call

summary(fit.gls.gj) # summary of the gls model; same as with IC


# without the grizzly and polar bear
summary(phylog.gls.fit(cbind(body.mass,hind.l.length), running.speed, Garland.Janis.Cov, exclude.tips=c("Tm","Ur"))) # as expected, not identically equal to IC
