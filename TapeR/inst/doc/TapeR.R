### R code from vignette source 'TapeR.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: TapeR.Rnw:48-49 (eval = FALSE)
###################################################
## install.packages("TapeR")


###################################################
### code chunk number 2: TapeR.Rnw:55-56
###################################################
library(TapeR)


###################################################
### code chunk number 3: TapeR.Rnw:60-61 (eval = FALSE)
###################################################
## ?TapeR


###################################################
### code chunk number 4: TapeR.Rnw:80-108
###################################################

#load example data
data(DxHx.df)

#prepare the data (could be defined in the function directly)
Id = DxHx.df[,"Id"]
x = DxHx.df[,"Hx"]/DxHx.df[,"Ht"]#calculate relative heights
y = DxHx.df[,"Dx"]

#plot the example data
plot(x,y,pch=as.character(Id),xlab="Relative height (m)", ylab="Diameter (cm)")


#define the relative knot positions and order of splines

knt_x = c(0.0, 0.1, 0.75, 1.0) # B-Spline knots: fix effects
ord_x = 4 # ord = order (4 = cubic)
knt_z = c(0.0, 0.1, 1.0); ord_z = 4 # B-Spline knots: rnd effects

#fit the model
taper.model <- TapeR_FIT_LME.f(Id, x, y, knt_x, ord_x, knt_z, ord_z,
IdKOVb = "pdSymm")

#save model parameters for documentation or dissimination
#parameters can be load()-ed and used to predict the taper
#or volume using one or several measured dbh
spruce.taper.pars <- taper.model$par.lme
#save(spruce.taper.pars, file="spruce.taper.pars.rdata")##uncomment to save


