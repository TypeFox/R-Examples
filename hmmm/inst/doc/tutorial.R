### R code from vignette source 'tutorial.rnw'

###################################################
### code chunk number 1: tutorial.rnw:181-182
###################################################
options(prompt = "R> ", continue = "+  ", width = 80, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: tutorial.rnw:184-187
###################################################
library("hmmm")
data("accident", package = "hmmm")
accident[1:20,]


###################################################
### code chunk number 3: tutorial.rnw:198-199
###################################################
y <- getnames(accident, st = 9)


###################################################
### code chunk number 4: tutorial.rnw:201-204
###################################################
count<-cbind(row.names(y)[1:3],y[1:3])
colnames(count)<-c("cell names", "counts")
print(count,quote=F)


###################################################
### code chunk number 5: tutorial.rnw:219-221
###################################################
margin <- marg.list(c("marg-marg-b-b", "b-marg-b-b", 
"marg-b-b-b", "b-b-b-b"))


###################################################
### code chunk number 6: tutorial.rnw:226-229
###################################################
model <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"))
model


###################################################
### code chunk number 7: tutorial.rnw:243-246
###################################################
modelB <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"),
sel = c(12:13, 14:17))


###################################################
### code chunk number 8: tutorial.rnw:250-252
###################################################
modB <- hmmm.mlfit(y, modelB)
modB


###################################################
### code chunk number 9: tutorial.rnw:256-257 (eval = FALSE)
###################################################
## print(modB, aname = "model B", printflag = TRUE)


###################################################
### code chunk number 10: tutorial.rnw:261-262 (eval = FALSE)
###################################################
## summary(modB)


###################################################
### code chunk number 11: tutorial.rnw:280-284
###################################################
modelA <- hmmm.model(marg = margin, lev = c(3, 4, 3, 2),
names = c("Type", "Time", "Age", "Hour"), sel = c(12:13, 14:17),
formula = ~ Type * Age * Hour + Time * Age * Hour + Type : Time)
modA <- hmmm.mlfit(y, modelA)


###################################################
### code chunk number 12: tutorial.rnw:288-289
###################################################
anova(modA, modB)


###################################################
### code chunk number 13: tutorial.rnw:296-300 (eval = FALSE)
###################################################
## modellog <- loglin.model(lev = c(3, 4, 3, 2),
## formula = ~ Type * Age * Hour + Time * Age * Hour + Type : Time,
## names = c("Type", "Time", "Age", "Hour"))
## modlog <- hmmm.mlfit(y, modellog)


###################################################
### code chunk number 14: tutorial.rnw:335-337
###################################################
margin <- marg.list(c("marg-marg-l-l", "g-marg-l-l", 
"marg-g-l-l", "g-g-l-l"))


###################################################
### code chunk number 15: tutorial.rnw:348-351
###################################################
model <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4), 
names = c("In", "Sa", "Co", "Ho"))
model


###################################################
### code chunk number 16: tutorial.rnw:358-364
###################################################
model1 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 26:27, 34:39))
data("madsen", package = "hmmm")
y <- getnames(madsen, st = 6)
mod1 <- hmmm.mlfit(y, model1)
mod1


###################################################
### code chunk number 17: tutorial.rnw:368-372
###################################################
model2 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 34:39))
mod2 <- hmmm.mlfit(y, model2)
mod2


###################################################
### code chunk number 18: tutorial.rnw:381-385
###################################################
model3 <- hmmm.model(marg = margin, lev = c(3, 3, 2, 4),
names = c("In", "Sa", "Co", "Ho"), sel = c(18:23, 34:39, 44:71))
mod3 <- hmmm.mlfit(y, model3)
mod3


###################################################
### code chunk number 19: tutorial.rnw:429-430
###################################################
marginals <- marg.list(c("r-marg", "marg-r", "r-r"))


###################################################
### code chunk number 20: tutorial.rnw:443-445
###################################################
rec1 <- matrix(c(-1, -1,  1,
                 -1,  1,  0), 2, 3, byrow = TRUE)


###################################################
### code chunk number 21: tutorial.rnw:454-460
###################################################
rec2 <- matrix(c(-1, -1, -1,  0,  1,  1,  1,
                 -1, -1, -1,  1,  0,  0,  0,
                  1, -1,  0,  0,  0,  0,  0,
                  0, -1,  1,  0,  0,  0,  0,
                  0,  0,  0,  0,  1, -1,  0,
                  0,  0,  0,  0,  0, -1,  1), 6, 7, byrow = TRUE)


###################################################
### code chunk number 22: tutorial.rnw:468-469
###################################################
rec <- recursive(rec1, rec2)


###################################################
### code chunk number 23: tutorial.rnw:477-480
###################################################
model <- hmmm.model(marg = marginals, lev = c(3, 7), 
names = c("Rel", "Pol"), cocacontr = rec)
model


###################################################
### code chunk number 24: tutorial.rnw:513-517
###################################################
Emat <- cbind(matrix(0, 2, 4), matrix(c(1, 0, 0, 1, 0, -1, -1, 0), 2, 4), 
matrix(0, 2, 12))
modelE <- hmmm.model(marg = marginals, lev = c(3, 7), 
names = c("Rel", "Pol"), cocacontr = rec, E = Emat)


###################################################
### code chunk number 25: tutorial.rnw:522-526
###################################################
data("relpol", package = "hmmm")
y <- getnames(relpol, st = 4)
modE <- hmmm.mlfit(y, modelE)
print(modE)


###################################################
### code chunk number 26: tutorial.rnw:542-544
###################################################
data("depression", package="hmmm")
y <- getnames(depression, st = 9)


###################################################
### code chunk number 27: tutorial.rnw:568-573
###################################################
margin <- marg.list(c("marg-marg-marg-b-b","b-marg-marg-b-b",
"marg-b-marg-b-b", "marg-marg-b-b-b","b-b-b-b-b"))
name <- c("R3","R2","R1","T","D")
modelsat<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name)
modelsat


###################################################
### code chunk number 28: tutorial.rnw:591-603
###################################################
A1<-matrix(c(
   0,0,0,1,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,1,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,1,
   0,1,0,0,0,0,0,0,0,-1,0,0,
   0,0,0,0,0,1,0,0,0,-1,0,0,
   0,0,1,0,0,0,0,0,0,0,-1,0,
   0,0,0,0,0,0,1,0,0,0,-1,0,
   1,0,0,0,-2,0,0,0,1,0,0,0
   ),8,12,byrow=TRUE)

E1<-cbind(matrix(0,8,3), A1, matrix(0,8,16))


###################################################
### code chunk number 29: tutorial.rnw:608-612
###################################################
model1<-hmmm.model(marg = margin, lev =c(2,2,2,2,2), names = name, E = E1)

fitmod1 <- hmmm.mlfit(y, model1)
fitmod1


###################################################
### code chunk number 30: tutorial.rnw:617-628
###################################################
A2<-matrix(c(
    0,0,0,1,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,1,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,1,
    0,1,0,0,0,-2,0,0,0,1,0,0,
    0,0,1,0,0,0,0,0,0,0,-1,0,
    0,0,0,0,0,0,1,0,0,0,-1,0,
    1,0,0,0,-2,0,0,0,1,0,0,0),
    7,12,byrow=TRUE)

E2<-cbind(matrix(0,7,3), A2, matrix(0,7,16))


###################################################
### code chunk number 31: tutorial.rnw:633-637
###################################################
model2<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name, E = E2)

fitmod2 = hmmm.mlfit(y, model2)
fitmod2


###################################################
### code chunk number 32: tutorial.rnw:647-652
###################################################
model3<-hmmm.model(marg = margin, lev = c(2,2,2,2,2), names = name, E = E2,
                   formula=~R1*R2*T*D+R3*R2*T*D )

fitmod3 = hmmm.mlfit(y, model3)
fitmod3


###################################################
### code chunk number 33: tutorial.rnw:678-679
###################################################
marginals <- marg.list(c("b-marg", "marg-g", "b-g"))


###################################################
### code chunk number 34: tutorial.rnw:693-698
###################################################
al <- list(
Type = ~ Type * (Age + Hour),
Time = ~ Time * (Age + Hour),
Type.Time = ~ Type.Time * (Age + Hour)
)


###################################################
### code chunk number 35: tutorial.rnw:713-716
###################################################
model <- hmmm.model.X(marg = marginals, lev = c(3, 4), 
names = c("Type", "Time"), Formula = al, strata = c(3, 2), 
fnames = c("Age", "Hour"))


###################################################
### code chunk number 36: tutorial.rnw:722-726
###################################################
data("accident", package = "hmmm")
y <- getnames(accident, st = 9)
mod1 <- hmmm.mlfit(y, model)
mod1


###################################################
### code chunk number 37: tutorial.rnw:732-733 (eval = FALSE)
###################################################
## summary(mod1)


###################################################
### code chunk number 38: tutorial.rnw:745-750
###################################################
alind <- list(
Type = ~ Type * Age + Type * Hour,
Time = ~ Time * Age + Time * Hour,
Type.Time = "zero"
)


###################################################
### code chunk number 39: tutorial.rnw:762-767
###################################################
alpar <- list(
Type = ~ Type + Age + Hour,
Time = ~ Time + Age + Hour,
Type.Time = ~ Type.Time + Age + Hour
)


###################################################
### code chunk number 40: tutorial.rnw:800-804
###################################################
data("polbirth", package = "hmmm")
y <- getnames(polbirth)
marginals <- marg.list(c("g-marg", "marg-l", "g-l"))
names <- c("Politics", "Birth")


###################################################
### code chunk number 41: tutorial.rnw:810-811
###################################################
ineq <- list(marg = c(1, 2), int = list(c(1, 2)), types = c("g", "l"))


###################################################
### code chunk number 42: tutorial.rnw:817-819
###################################################
model <- hmmm.model(marg = marginals, dismarg = ineq, lev = c(7, 4), 
names = names)


###################################################
### code chunk number 43: tutorial.rnw:829-830
###################################################
mlr <- hmmm.mlfit(y, model, noineq = FALSE)


###################################################
### code chunk number 44: tutorial.rnw:836-837
###################################################
msat <- hmmm.mlfit(y, model)


###################################################
### code chunk number 45: tutorial.rnw:843-846
###################################################
model0 <- hmmm.model(marg = marginals, lev = c(7, 4), sel = c(10:27), 
names = names)
mnull <- hmmm.mlfit(y, model0)


###################################################
### code chunk number 46: tutorial.rnw:856-857
###################################################
test <- hmmm.chibar(nullfit = mnull, disfit = mlr, satfit = msat)


###################################################
### code chunk number 47: tutorial.rnw:876-877
###################################################
test


###################################################
### code chunk number 48: tutorial.rnw:898-899
###################################################
y <- matrix(c(104, 24, 65, 76, 146, 30, 50, 9, 166), 9, 1)


###################################################
### code chunk number 49: tutorial.rnw:910-912
###################################################
Zmat <- kronecker(diag(3), matrix(1, 3, 1))
ZFmat <- kronecker(diag(3), matrix(1, 3, 1))[,3]


###################################################
### code chunk number 50: tutorial.rnw:924-931
###################################################
Gini <- function(m) {
A<-matrix(m,3,3,byrow=TRUE)
 GNum<-rowSums(A^2)
 GDen<-rowSums(A)^2
 G<-GNum/GDen
 c(G[1], G[3]) - c(G[2], G[1])
 }


###################################################
### code chunk number 51: tutorial.rnw:941-942
###################################################
mod_eq <- mphineq.fit(y, Z = Zmat, ZF = ZFmat, h.fct = Gini)


###################################################
### code chunk number 52: tutorial.rnw:951-952
###################################################
mod_ineq <- mphineq.fit(y, Z = Zmat, ZF = ZFmat, d.fct = Gini)


###################################################
### code chunk number 53: tutorial.rnw:958-960
###################################################
mod_sat <- mphineq.fit(y, Z = Zmat, ZF = ZFmat)
hmmm.chibar(nullfit = mod_eq, disfit = mod_ineq, satfit = mod_sat)


