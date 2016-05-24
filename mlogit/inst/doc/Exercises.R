### R code from vignette source 'Exercises.rnw'

###################################################
### code chunk number 1: Exercises.rnw:48-49
###################################################
options(prompt= "R> ", useFancyQuotes = FALSE)


###################################################
### code chunk number 2: Exercises.rnw:124-130
###################################################
library("mlogit")
#source("~/Forge/mlogit/chargement.R")
data("Heating", package = "mlogit")
H <- mlogit.data(Heating, shape="wide", choice="depvar", varying=c(3:12))
m <- mlogit(depvar~ic+oc|0, H)
summary(m)


###################################################
### code chunk number 3: Exercises.rnw:147-148
###################################################
apply(fitted(m, outcome=FALSE), 2, mean)


###################################################
### code chunk number 4: Exercises.rnw:170-171
###################################################
coef(m)["oc"]/coef(m)["ic"]


###################################################
### code chunk number 5: Exercises.rnw:218-222
###################################################
H$lcc=H$ic+H$oc/0.12
mlcc <- mlogit(depvar~lcc|0, H)
lrtest(m, mlcc)
qchisq(0.05, df = 1, lower.tail = FALSE)


###################################################
### code chunk number 6: Exercises.rnw:244-247
###################################################
mc <- mlogit(depvar~ic+oc, H, reflevel = 'hp')
summary(mc)
apply(fitted(mc, outcome = FALSE), 2, mean)


###################################################
### code chunk number 7: Exercises.rnw:257-261
###################################################
wtp <- coef(mc)["oc"]/coef(mc)["ic"]
wtp
r <- 1/wtp
r


###################################################
### code chunk number 8: Exercises.rnw:296-297
###################################################
update(mc, reflevel="gr")


###################################################
### code chunk number 9: Exercises.rnw:310-312
###################################################
mi <- mlogit(depvar~oc+I(ic/income), H)
summary(mi)


###################################################
### code chunk number 10: Exercises.rnw:324-325
###################################################
mi2 <- mlogit(depvar~oc+ic|income, H, reflevel="hp")


###################################################
### code chunk number 11: Exercises.rnw:338-341
###################################################
lrtest(mc, mi2)
waldtest(mc, mi2)
scoretest(mc, mi2)


###################################################
### code chunk number 12: Exercises.rnw:358-367
###################################################
X <- model.matrix(mc)
alt <- index(H)$alt
chid <- index(H)$chid
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- tapply(eXb, chid, sum)
P <- eXb / SeXb[chid]
P <- matrix(P, ncol = 5, byrow = TRUE)
head(P)
apply(P, 2, mean)


###################################################
### code chunk number 13: Exercises.rnw:377-378
###################################################
apply(fitted(mc, outcome = FALSE), 2, mean)


###################################################
### code chunk number 14: Exercises.rnw:391-394
###################################################
Hn <- H
Hn[Hn$alt == "hp", "ic"] <- 0.9 * Hn[Hn$alt == "hp", "ic"]
apply(predict(mc, newdata = Hn), 2, mean)


###################################################
### code chunk number 15: Exercises.rnw:417-431
###################################################
X <- model.matrix(mc)
Xn <- X[alt == "ec",]
Xn[, "ic"] <- Xn[, "ic"] + 200
Xn[, "oc"] <- Xn[, "oc"] * 0.75
unchid <- unique(index(H)$chid)
rownames(Xn) <- paste(unchid, 'new', sep = ".")
chidb <- c(chid, unchid)
X <- rbind(X, Xn)
X <- X[order(chidb), ]
eXb <- as.numeric(exp(X %*% coef(mc)))
SeXb <- as.numeric(tapply(eXb, sort(chidb), sum))
P <- eXb / SeXb[sort(chidb)]
P <- matrix(P, ncol = 6, byrow = TRUE)
apply(P, 2, mean)


###################################################
### code chunk number 16: Exercises.rnw:497-517
###################################################
library("mlogit")
data("HC", package = "mlogit")
HC <- mlogit.data(HC, varying = c(2:8, 10:16), choice = "depvar", shape = "wide")
cooling.modes <- index(HC)$alt %in% c('gcc', 'ecc', 'erc', 'hpc')
room.modes <- index(HC)$alt %in% c('erc', 'er')
# installation / operating costs for cooling are constants, 
# only relevant for mixed systems
HC$icca[!cooling.modes] <- 0
HC$occa[!cooling.modes] <- 0
# create income variables for two sets cooling and rooms
HC$inc.cooling <- HC$inc.room <- 0
HC$inc.cooling[cooling.modes] <- HC$income[cooling.modes]
HC$inc.room[room.modes] <- HC$income[room.modes]
# create an intercet for cooling modes
HC$int.cooling <- as.numeric(cooling.modes)
# estimate the model with only one nest elasticity
nl <- mlogit(depvar ~ ich + och +icca + occa + inc.room + inc.cooling + int.cooling | 0, HC,
             nests = list(cooling = c('gcc','ecc','erc','hpc'), 
             other = c('gc', 'ec', 'er')), un.nest.el = TRUE)
summary(nl)


###################################################
### code chunk number 17: Exercises.rnw:537-538
###################################################
 (coef(nl)['iv']-1)/sqrt(vcov(nl)['iv', 'iv'])


###################################################
### code chunk number 18: Exercises.rnw:545-548
###################################################
# First estimate the multinomial logit model
ml <- update(nl, nests = NULL)
lrtest(nl, ml)


###################################################
### code chunk number 19: Exercises.rnw:559-562
###################################################
nl2 <- update(nl, nests = list(central = c('ec', 'ecc', 'gc', 'gcc', 'hpc'), 
                    room = c('er', 'erc')))
summary(nl2)


###################################################
### code chunk number 20: Exercises.rnw:579-581
###################################################
 (coef(nl2)['iv']-1)/sqrt(vcov(nl2)['iv', 'iv'])
lrtest(nl2, ml)


###################################################
### code chunk number 21: Exercises.rnw:591-593
###################################################
logLik(nl)
logLik(nl2)


###################################################
### code chunk number 22: Exercises.rnw:602-603
###################################################
nl3 <- update(nl, un.nest.el = FALSE)


###################################################
### code chunk number 23: Exercises.rnw:624-625
###################################################
lrtest(nl, nl3)


###################################################
### code chunk number 24: Exercises.rnw:643-646
###################################################
nl4 <- update(nl, nests=list(n1 = c('gcc', 'ecc', 'erc'), n2 = c('hpc'),
                    n3 = c('gc', 'ec', 'er')))
summary(nl4)


###################################################
### code chunk number 25: Exercises.rnw:694-698
###################################################
library("mlogit")
data("Electricity", package = "mlogit")
Electr <- mlogit.data(Electricity, id="id", choice="choice", 
                      varying=3:26, shape="wide", sep="")


###################################################
### code chunk number 26: Exercises.rnw:700-703 (eval = FALSE)
###################################################
## Elec.mxl <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
##               rpar=c(pf='n', cl='n', loc='n', wk='n', tod='n', seas='n'), 
##               R=100, halton=NA, print.level=0, panel=TRUE)


###################################################
### code chunk number 27: Exercises.rnw:705-706
###################################################
data("Examples", package = "mlogit")


###################################################
### code chunk number 28: Exercises.rnw:708-709
###################################################
summary(Elec.mxl)


###################################################
### code chunk number 29: Exercises.rnw:716-717
###################################################
coef(Elec.mxl)['cl']/coef(Elec.mxl)['pf']


###################################################
### code chunk number 30: Exercises.rnw:730-731
###################################################
pnorm(-coef(Elec.mxl)['cl']/coef(Elec.mxl)['sd.cl'])


###################################################
### code chunk number 31: Exercises.rnw:759-760
###################################################
pnorm(-coef(Elec.mxl)['pf']/coef(Elec.mxl)['sd.pf'])


###################################################
### code chunk number 32: Exercises.rnw:768-771 (eval = FALSE)
###################################################
## Elec.mxl2 <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
##                    rpar=c(cl='n', loc='n', wk='n', tod='n', seas='n'), 
##                    R=100, halton=NA, print.level=0, panel=TRUE)


###################################################
### code chunk number 33: Exercises.rnw:773-774
###################################################
summary(Elec.mxl2)


###################################################
### code chunk number 34: Exercises.rnw:785-786 (eval = FALSE)
###################################################
## Elec.mxl3 <- update(Elec.mxl, rpar=c(cl='n', loc='n', wk='u', tod='n', seas='n'))


###################################################
### code chunk number 35: Exercises.rnw:790-793
###################################################
summary(Elec.mxl3)
rpar(Elec.mxl3, 'wk')
summary(rpar(Elec.mxl3, 'wk'))


###################################################
### code chunk number 36: Exercises.rnw:795-796
###################################################
plot(rpar(Elec.mxl3, 'wk'))


###################################################
### code chunk number 37: Exercises.rnw:835-841 (eval = FALSE)
###################################################
## Electr <- mlogit.data(Electricity, id="id", choice="choice", 
##                       varying=3:26, shape="wide", sep="",
##                       opposite=c('tod', 'seas'))
## Elec.mxl4 <- mlogit(choice~pf+cl+loc+wk+tod+seas|0, Electr, 
##               rpar=c(cl='n', loc='n', wk='u', tod='ln', seas='ln'), 
##               R=100, halton=NA, print.level=0, panel=TRUE)


###################################################
### code chunk number 38: Exercises.rnw:843-844
###################################################
summary(Elec.mxl4)


###################################################
### code chunk number 39: Exercises.rnw:846-847
###################################################
plot(rpar(Elec.mxl4, 'seas'))


###################################################
### code chunk number 40: Exercises.rnw:855-856 (eval = FALSE)
###################################################
## Elec.mxl5 <- update(Elec.mxl4, correlation = TRUE)


###################################################
### code chunk number 41: Exercises.rnw:858-863
###################################################
summary(Elec.mxl5)
cor.mlogit(Elec.mxl5)
lrtest(Elec.mxl5, Elec.mxl4)
waldtest(Elec.mxl5, correlation = FALSE)
scoretest(Elec.mxl4, correlation = TRUE)


###################################################
### code chunk number 42: Exercises.rnw:918-920
###################################################
data("Mode", package="mlogit")
Mo <- mlogit.data(Mode, choice='choice', shape='wide', varying=c(2:9))


###################################################
### code chunk number 43: Exercises.rnw:923-924 (eval = FALSE)
###################################################
## p1 <- mlogit(choice~cost+time, Mo, seed = 20, R = 100, probit = TRUE)


###################################################
### code chunk number 44: Exercises.rnw:927-928
###################################################
summary(p1)


###################################################
### code chunk number 45: Exercises.rnw:932-934
###################################################
L1 <- matrix(0, 3, 3)
L1[!upper.tri(L1)] <- c(1, coef(p1)[6:10])


###################################################
### code chunk number 46: Exercises.rnw:939-940
###################################################
L1 %*% t(L1)


###################################################
### code chunk number 47: Exercises.rnw:977-978 (eval = FALSE)
###################################################
## p2 <- mlogit(choice~cost+time, Mo, seed = 21, R = 100, probit = TRUE)


###################################################
### code chunk number 48: Exercises.rnw:981-982
###################################################
coef(p2)


###################################################
### code chunk number 49: Exercises.rnw:1001-1002
###################################################
actShares <- with(Mo, tapply(choice, alt, mean))


###################################################
### code chunk number 50: Exercises.rnw:1008-1011
###################################################
predShares <- apply(fitted(p1, outcome = FALSE), 2, mean)
predShares
sum(predShares)


###################################################
### code chunk number 51: Exercises.rnw:1032-1036
###################################################
Mo2 <- Mo
Mo2[Mo2$alt == 'car', 'cost'] <- Mo2[Mo2$alt == 'car', 'cost'] * 2
newShares <- apply(predict(p1, newdata = Mo2), 2, mean)
cbind(original = actShares, new = newShares, change = round((newShares - actShares) / actShares * 100))


