###################################################
### chunk number 1: 
###################################################
options("digits"=4)
options("width"=80)


###################################################
### chunk number 2: 
###################################################
require(cem)
data(LeLonde)


###################################################
### chunk number 3: 
###################################################
Le <- data.frame(na.omit(LeLonde))


###################################################
### chunk number 4: 
###################################################
tr <- which(Le$treated==1)
ct <- which(Le$treated==0)
ntr <- length(tr)
nct <- length(ct)


###################################################
### chunk number 5: 
###################################################
mean(Le$re78[tr]) - mean(Le$re78[ct])


###################################################
### chunk number 6: 
###################################################
vars <- c("age", "education", "black", "married", "nodegree", "re74",
"re75", "hispanic", "u74", "u75","q1")


###################################################
### chunk number 7: 
###################################################
L1 <- L1.meas(Le$treated, Le[vars])$L1 


###################################################
### chunk number 8: 
###################################################
imbalance(group=Le$treated, data=Le[vars])


###################################################
### chunk number 9: 
###################################################
mat <- cem(treatment = "treated", data = Le, drop = "re78")


###################################################
### chunk number 10: 
###################################################
mat


###################################################
### chunk number 11: 
###################################################
levels(Le$q1)


###################################################
### chunk number 12: 
###################################################
q1.grp <- list(c("strongly agree", "agree"), c("neutral","no opinion"), c("strongly disagree","disagree"))


###################################################
### chunk number 13: 
###################################################
table(Le$education)


###################################################
### chunk number 14: 
###################################################
educut <- c(0, 6.5, 8.5, 12.5, 17)


###################################################
### chunk number 15: 
###################################################
mat1 <- cem(treatment = "treated", data = Le, drop = "re78", 
cutpoints = list(education=educut), grouping=list(q1=q1.grp))
mat1


###################################################
### chunk number 16: 
###################################################
mat$breaks$education


###################################################
### chunk number 17: 
###################################################
mat1$breaks$education


###################################################
### chunk number 18: 
###################################################
cem("treated", Le, cutpoints = list(age=10), drop="re78", grouping=list(q1=q1.grp))
cem("treated", Le, cutpoints = list(age=6), drop="re78", grouping=list(q1=q1.grp))
cem("treated", Le, cutpoints = list(age=3), drop="re78", grouping=list(q1=q1.grp))


###################################################
### chunk number 19: 
###################################################
tab <- relax.cem(mat, Le, depth=1, perc=0.3) 


###################################################
### chunk number 20: 
###################################################
pdf("coarsen1.pdf", width=9, height=6, pointsize=10)
plot(tab,perc=0.3)
invisible(dev.off())


###################################################
### chunk number 21: 
###################################################
plot(tab, group="1", perc=0.35,unique=TRUE)


###################################################
### chunk number 22: 
###################################################
pdf("coarsen2.pdf", width=9, height=6, pointsize=10)
plot(tab, group="1", perc=0.35,unique=TRUE)
invisible(dev.off())


###################################################
### chunk number 23: 
###################################################
mat <- cem(treatment="treated",data=Le, drop="re78")
mat
mat$k2k


###################################################
### chunk number 24: 
###################################################
mat2 <- k2k(mat, Le, "euclidean", 1)
mat2
mat2$k2k


###################################################
### chunk number 25: 
###################################################
data(LL)
mat <- cem(treatment="treated", data=LL, drop="re78")
est <- att(mat, re78 ~ treated, data = LL)
est


###################################################
### chunk number 26: 
###################################################
est2 <- att(mat, re78 ~ treated + re74, data = LL)
est2


###################################################
### chunk number 27: 
###################################################
att(mat, re78 ~ treated + re74 , data = LL, model="linear")


###################################################
### chunk number 28: 
###################################################
att(mat, re78 ~ treated + re74 , data = LL, model="linear-RE")


###################################################
### chunk number 29: 
###################################################
att(mat, re78 ~ treated + re74 , data = LL, model="forest")


###################################################
### chunk number 30: 
###################################################
att(mat, re78 ~ treated + re74 , data = LL, model="linear", extra=TRUE)
att(mat, re78 ~ treated + re74 , data = LL, model="linear-RE", extra=TRUE)
att(mat, re78 ~ treated + re74 , data = LL, model="rf", extra=TRUE)


###################################################
### chunk number 31: 
###################################################
est3 <- att(mat, re78 ~ treated + re74 , data = LL)
est3
plot(est3, mat, LL, vars=c("education", "age", "re74", "re75"))


###################################################
### chunk number 32: 
###################################################
pdf("teff.pdf", width=9, height=6, pointsize=10)
est3 <- att(mat, re78 ~ treated + re74 + re75, data = LL)
plot(est3, mat, LL, vars=c("education", "age", "re74", "re75"))
invisible(dev.off())


###################################################
### chunk number 33: 
###################################################
mat3 <- cem("treated", LeLonde, drop="re78", cutpoints = mat$breaks, grouping=list(q1=q1.grp))
mat3


###################################################
### chunk number 34: 
###################################################
mat4 <- cem("treated", Le, drop="re78", cutpoints = mat$breaks, grouping=list(q1=q1.grp))
mat4


###################################################
### chunk number 35: 
###################################################
summary(LeLonde)


###################################################
### chunk number 36: 
###################################################
require(Amelia)
set.seed(123)
imputed <- amelia(LeLonde,noms=c("black","hispanic","treated","married","nodegree",
"u74","u75","q1"))
imputed <- imputed$imputations[1:5]


###################################################
### chunk number 37: 
###################################################
mat2 <- cem("treated", datalist=imputed, drop="re78", data=LeLonde, grouping=list(q1=q1.grp))
mat2


###################################################
### chunk number 38: 
###################################################
out <- att(mat2, re78 ~ treated, data=imputed)
out


###################################################
### chunk number 39: 
###################################################
data(LL)

# cem match: automatic bin choice
mat <- cem(data=LL, drop="re78")

# we want a set of paired units
psample <- pair(mat, data=LL)


###################################################
### chunk number 40: 
###################################################
table(psample$paired)


###################################################
### chunk number 41: 
###################################################
psample$paired[1:100]


###################################################
### chunk number 42: 
###################################################
table(psample$full.paired)
psample$full.paired[1:10]


###################################################
### chunk number 43: 
###################################################
# cem match: automatic bin choice, we drop one row from the data set
mat1 <- cem(data=LL[-1,], drop="re78")

# we want a set of paired units but we have an odd number of units in the data
psample <- pair(mat1, data=LL[-1,])
table(psample$full.paired)


