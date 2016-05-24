### R code from vignette source 'gplm-examples.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gplm-examples.Rnw:818-820
###################################################
data(airquality)
attach(airquality)


###################################################
### code chunk number 2: gplm-examples.Rnw:824-827
###################################################
library(gplm)
fh <- kde(Wind)
plot(fh, type="l", main="Kernel density estimate (KDE)")


###################################################
### code chunk number 3: gplm-examples.Rnw:833-834
###################################################
plot(fh, type="l", main="Kernel density estimate (KDE)")


###################################################
### code chunk number 4: gplm-examples.Rnw:844-845
###################################################
fh$bandwidth


###################################################
### code chunk number 5: gplm-examples.Rnw:851-854
###################################################
fh.10 <- kde(Wind, grid=c(10,15))
points(fh.10, col="red", pch=19)
title("KDE with estimates at x=10, 15 (in red)")


###################################################
### code chunk number 6: gplm-examples.Rnw:861-864
###################################################
plot(fh, type="l")
points(fh.10, col="red", pch=19)
title("KDE with estimates at x=10, 15 (in red)")


###################################################
### code chunk number 7: gplm-examples.Rnw:871-873
###################################################
fh <- kde(Wind, bandwidth=3, kernel="epanechnikov")
fh$bandwidth


###################################################
### code chunk number 8: gplm-examples.Rnw:879-881
###################################################
plot(fh, type="l")
title("KDE with uniform kernel and bandwidth=3")


###################################################
### code chunk number 9: gplm-examples.Rnw:887-888
###################################################
fh.biv <- kde(cbind(Wind,Temp))


###################################################
### code chunk number 10: gplm-examples.Rnw:891-897
###################################################
Wind.grid <- unique(fh.biv$x[,1])  ## extract grids
Temp.grid <- unique(fh.biv$x[,2])  
o <- order(fh.biv$x[,2],fh.biv$x[,1])  ## order by 2nd column
fh2 <- matrix(fh.biv$y[o],length(Wind.grid),length(Temp.grid))
persp(Wind.grid,Temp.grid,fh2,main="bivariate KDE",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 11: gplm-examples.Rnw:903-905
###################################################
persp(Wind.grid,Temp.grid,fh2,main="bivariate KDE",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 12: gplm-examples.Rnw:910-917
###################################################
Wind.grid <- seq(min(Wind),max(Wind),length=20)  ## define grid
Temp.grid <- seq(min(Temp),max(Temp),length=40)
fh.biv <- kde(cbind(Wind,Temp), grid=create.grid(list(Wind.grid,Temp.grid)))
o <- order(fh.biv$x[,2],fh.biv$x[,1])  ## order by 2nd column
fh2a <- matrix(fh.biv$y[o],length(Wind.grid),length(Temp.grid))
persp(Wind.grid,Temp.grid,fh2a,main="bivariate KDE",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 13: gplm-examples.Rnw:923-925
###################################################
persp(Wind.grid,Temp.grid,fh2a,main="bivariate KDE, different grids",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 14: gplm-examples.Rnw:930-931
###################################################
contour(Wind.grid,Temp.grid,fh2a, main="KDE Contours")


###################################################
### code chunk number 15: gplm-examples.Rnw:935-936
###################################################
contour(Wind.grid,Temp.grid,fh2a, main="KDE Contours")


###################################################
### code chunk number 16: gplm-examples.Rnw:943-947
###################################################
mh <- kreg(Wind, Temp)
plot(Wind,Temp, col="grey", pch="+")
lines(mh, col="blue", lwd=2)
title("Data and Nadaraya--Watson regression")


###################################################
### code chunk number 17: gplm-examples.Rnw:954-957
###################################################
plot(Wind,Temp, col="grey", pch="+")
lines(mh, col="blue", lwd=2)
title("Data and Nadaraya--Watson regression")


###################################################
### code chunk number 18: gplm-examples.Rnw:963-974
###################################################
airquality2 <- airquality[!is.na(Ozone),] ## delete NA's
detach(airquality)  ## detach previous data
attach(airquality2)
bandwidth <- bandwidth.scott(cbind(Wind,Temp))
mh.biv <- kreg(cbind(Wind,Temp),Ozone, bandwidth=bandwidth)
Wind.grid <- unique(mh.biv$x[,1])  ## extract grids
Temp.grid <- unique(mh.biv$x[,2])  
o <- order(mh.biv$x[,2],mh.biv$x[,1])  ## order by 2nd column
mh2 <- matrix(mh.biv$y[o],length(Wind.grid),length(Temp.grid))
persp(Wind.grid,Temp.grid,mh2,main="bivariate KDE",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 19: gplm-examples.Rnw:980-982
###################################################
persp(Wind.grid,Temp.grid,mh2,main="bivariate KDE",
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)


###################################################
### code chunk number 20: gplm-examples.Rnw:992-996
###################################################
library(AER)
data(CPS1985)
str(CPS1985)  ## show data structure
attach(CPS1985)


###################################################
### code chunk number 21: gplm-examples.Rnw:1002-1008
###################################################
library(gplm)
bandwidth <- bandwidth.scott(experience)
gh <- kgplm(x=cbind(gender,education),t=experience,y=log(wage),h=bandwidth)
o <- order(experience)  ## sort curve estimate by experience
plot(experience[o], gh$m[o], type="l")
title("PLM component function for experience")


###################################################
### code chunk number 22: gplm-examples.Rnw:1015-1017
###################################################
plot(experience[o], gh$m[o], type="l")
title("PLM component function for experience")


###################################################
### code chunk number 23: gplm-examples.Rnw:1023-1028
###################################################
exp.grid <- seq(min(experience),max(experience),length=200)
gh2 <- kgplm(x=cbind(gender,education),t=experience,y=log(wage),
             h=bandwidth,grid=exp.grid)
plot(exp.grid, gh2$m.grid, type="l", col="blue")
title("PLM component function for experience (on grid)")


###################################################
### code chunk number 24: gplm-examples.Rnw:1033-1035
###################################################
plot(exp.grid, gh2$m.grid, type="l", col="blue")
title("PLM component function for experience (on grid)")


###################################################
### code chunk number 25: gplm-examples.Rnw:1044-1048
###################################################
gs <- sgplm1(x=cbind(gender,education),t=experience,y=log(wage),df=8)
o <- order(experience)  ## sort curve estimate by experience
plot(experience[o], gs$m[o], type="l")
title("PLM component function for experience (sgplm1)")


###################################################
### code chunk number 26: gplm-examples.Rnw:1054-1056
###################################################
plot(experience[o], gs$m[o], type="l")
title("PLM component function for experience (sgplm1)")


###################################################
### code chunk number 27: gplm-examples.Rnw:1066-1079
###################################################
bandwidth <- 1.5*bandwidth.scott(cbind(education,experience))
edu.grid <- seq(min(education),max(education),length=25)
exp.grid <- seq(min(experience),max(experience),length=25)
grid  <- create.grid(list(edu.grid,exp.grid))

gh <- kgplm(x=(gender=="female"),t=cbind(education,experience),y=log(wage),
            h=bandwidth,grid=grid)

o <- order(grid[,2],grid[,1])
mh <- matrix(gh$m.grid[o],length(edu.grid),length(exp.grid))
persp(edu.grid,exp.grid,mh,
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
title("bivariate PLM component function")


###################################################
### code chunk number 28: gplm-examples.Rnw:1084-1087
###################################################
persp(edu.grid,exp.grid,mh,
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
title("bivariate PLM component function")


###################################################
### code chunk number 29: gplm-examples.Rnw:1096-1097
###################################################
detach(CPS1985)  ## detach previous data 


###################################################
### code chunk number 30: gplm-examples.Rnw:1103-1108
###################################################
library(AER)
data(Affairs)
str(Affairs)  ## show data structure
attach(Affairs)
y <- (affairs > 0)


###################################################
### code chunk number 31: gplm-examples.Rnw:1115-1122
###################################################
library(gplm)
bandwidth <- 1.5*bandwidth.scott(age)
age.grid <- seq(min(age),max(age),length=200)
gh <- kgplm(x=cbind(gender,education,yearsmarried),t=age,y=y,h=bandwidth,
            grid=age.grid,family="bernoulli",link="logit")
plot(age.grid, gh$m.grid, type="l")
title("GPLM component function for age")


###################################################
### code chunk number 32: gplm-examples.Rnw:1127-1129
###################################################
plot(age.grid, gh$m.grid, type="l")
title("GPLM component function for age")


###################################################
### code chunk number 33: gplm-examples.Rnw:1134-1141
###################################################
gs <- sgplm1(x=cbind(gender,education,yearsmarried),t=age,y=y,df=7,
             grid=age.grid,family="bernoulli",link="logit")
ylim <- range(gh$m.grid, gs$m.grid)
plot(age.grid, gh$m.grid, type="l", ylim=ylim)
lines(age.grid, gs$m.grid, col="seagreen")
title("GPLM component function for age (kgplm vs. sgplm1)")
legend("topright",c("kgplm","sgplm1"),lwd=1,col=c("black","seagreen"))


###################################################
### code chunk number 34: gplm-examples.Rnw:1146-1151
###################################################
ylim <- range(gh$m.grid, gs$m.grid)
plot(age.grid, gh$m.grid, type="l", ylim=ylim)
lines(age.grid, gs$m.grid, col="seagreen")
title("GPLM component function for age (kgplm vs. sgplm1)")
legend("topright",c("kgplm","sgplm1"),lwd=1,col=c("black","seagreen"))


###################################################
### code chunk number 35: biv.gplm
###################################################
bandwidth <- 1.5*bandwidth.scott(cbind(age,yearsmarried))
age.grid <- seq(min(age),max(age),length=25)
ym.grid <- seq(min(yearsmarried),max(yearsmarried),length=25)
grid  <- create.grid(list(age.grid,ym.grid))

gh <- kgplm(x=(gender=="female"),t=cbind(age,yearsmarried),y=y,
            h=bandwidth,grid=grid,family="bernoulli",link="logit")

o <- order(grid[,2],grid[,1])
mh <- matrix(gh$m.grid[o],length(age.grid),length(ym.grid))
persp(age.grid,ym.grid,mh,
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
title("bivariate GPLM component function")


###################################################
### code chunk number 36: gplm-examples.Rnw:1178-1181
###################################################
persp(age.grid,ym.grid,mh,
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
title("bivariate GPLM component function")


###################################################
### code chunk number 37: gplm-examples.Rnw:1187-1188
###################################################
bandwidth <- 1.5*bandwidth.scott(cbind(age,yearsmarried))
age.grid <- seq(min(age),max(age),length=25)
ym.grid <- seq(min(yearsmarried),max(yearsmarried),length=25)
grid  <- create.grid(list(age.grid,ym.grid))

gh <- kgplm(x=(gender=="female"),t=cbind(age,yearsmarried),y=y,
            h=bandwidth,grid=grid,family="bernoulli",link="logit")

o <- order(grid[,2],grid[,1])
mh <- matrix(gh$m.grid[o],length(age.grid),length(ym.grid))
persp(age.grid,ym.grid,mh,
      theta=30,phi=30,expand=0.5,col="lightblue",shade=0.5)
title("bivariate GPLM component function")


