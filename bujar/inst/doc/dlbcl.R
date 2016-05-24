### R code from vignette source 'dlbcl.Rnw'

###################################################
### code chunk number 1: dlbcl.Rnw:45-46
###################################################
options(prompt = "R> ", continue = " ", width = 90, digits =4, useFancyQuotes = FALSE)


###################################################
### code chunk number 2: dlbcl.Rnw:48-53
###################################################
library("bujar")
data("chop")
data("rchop")
###censoring rate in the CHOP data
sum(chop$status==0)/nrow(chop)


###################################################
### code chunk number 3: dlbcl.Rnw:55-57
###################################################
rchop <- subset(rchop, select=colnames(rchop)%in% colnames(chop))
chop$survtime <- chop$survtime + 1 ### add 1 for log-transformation


###################################################
### code chunk number 4: lin-fit (eval = FALSE)
###################################################
## set.seed(123)
## res.lin <- bujar(y=log(chop[,1]), cens=chop[,2], x=chop[,-(1:2)], tuning=TRUE, 
## cv=TRUE, mstop=1000)
## ###number of genes selected with BJ-LS
## sum(res.lin$xselect==1)
## coef.bj <- coef(res.lin)
## ###estimated non-zero coefficients (only list 10) 
## coef.bj[abs(coef.bj)>0][1:10]


###################################################
### code chunk number 5: dlbcl.Rnw:71-73
###################################################
library("survival")
cutyear <- 3


###################################################
### code chunk number 6: dlbcl.Rnw:75-81 (eval = FALSE)
###################################################
## pred.bj <- predict(res.lin, newx=rchop[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 7: dlbcl.Rnw:86-90 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 8: lin-twin-fit (eval = FALSE)
###################################################
## res.lin2 <- bujar(y=log(chop[,1]), cens=chop[,2], x=chop[,-(1:2)], tuning=TRUE, 
## cv=FALSE, mstop=1000, twin=TRUE, mstop2=100)
## ### number of genes selected with BJ-LS
## sum(res.lin2$xselect==1)
## coef.bj <- coef(res.lin2)
## coef.bj[abs(coef.bj)>0]
## pred.bj <- predict(res.lin2, newx=rchop[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 9: dlbcl.Rnw:115-119 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 10: preprocess (eval = FALSE)
###################################################
## library("rms")
## res <- rep(NA,ncol(chop))
## for(i in 3:ncol(chop)){
##   bjres <- try(bj(Surv(survtime, status) ~ chop[,i],data=chop, link="log"))
##   ###if BJ convergence fails, still included for further analysis
##   if(inherits(bjres, "try-error")) res[i] <- 1e-5 
##   else res[i] <- anova(bjres)[1,3]  #p-value
## }
## nsel <- 100
## ### select top nsel=100 genes with most significant p-values
## chop2 <- chop[, c(1, 2, sort.list(res,decreasing=FALSE)[1:nsel])]  
## rchop2 <- rchop[, c(1, 2, sort.list(res,decreasing=FALSE)[1:nsel])]
## colnames(chop2)[-(1:2)] <- colnames(rchop2)[-(1:2)]  <- 
## paste("x",colnames(chop2)[-(1:2)],sep="")
## detach(package:rms)


###################################################
### code chunk number 11: lasso-fit (eval = FALSE)
###################################################
## res.lasso <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="enet2", tuning=FALSE, whichlambda=20)
## ### how many genes selected by BJ-LASSO
## sum(res.lasso$xselect==1)
## ###estimated non-zero coefficients (only list 10) 
## coef.bj <- coef(res.lasso)
## coef.bj[abs(coef.bj)>0][1:10]
## pred.bj <- predict(res.lasso, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 12: dlbcl.Rnw:166-170 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 13: scad-fit (eval = FALSE)
###################################################
## res.scad <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)],  
## learner="snet", tuning=FALSE, whichlambda=20)
## ### how many genes selected by BJ-SCAD
## sum(res.scad$xselect==1)
## ###estimated non-zero coefficients (only list 10) 
## coef.bj <- coef(res.scad)
## coef.bj[abs(coef.bj)>0][1:10]
## pred.bj <- predict(res.scad, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 14: dlbcl.Rnw:196-200 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 15: sm-fit (eval = FALSE)
###################################################
## set.seed(123)
## res.ss <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="pspline", tuning=FALSE, cv=FALSE, mstop=100) 
## ### how many genes selected by BJ smoothing splines, only list 10
## sum(res.ss$xselect==1)
## colnames(res.ss$x)[res.ss$xselect==1][1:10]
## pred.bj <- predict(res.ss, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 16: dlbcl.Rnw:224-228 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 17: sm-twin-fit (eval = FALSE)
###################################################
## set.seed(123)
## res.ss2 <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="pspline", tuning=TRUE, cv=TRUE, mstop=100, twin=TRUE, mstop2=200)
## ### how many genes selected by BJ twin smoothing splines, only list 10
## sum(res.ss2$xselect==1) 
## colnames(res.ss2$x)[res.ss2$xselect==1][1:10]
## pred.bj <- predict(res.ss2, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 18: dlbcl.Rnw:251-255 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 19: tree-fit (eval = FALSE)
###################################################
## res.tree1 <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="tree",tuning=TRUE, cv=TRUE, mstop=1000, n.cores=2, rng=123)
## ###Number of genes selected with tree, only list 10
## sum(res.tree1$xselect==1)
## colnames(res.tree1$x)[res.tree1$xselect==1][1:10]
## pred.bj <- predict(res.tree1, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 20: dlbcl.Rnw:278-282 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 21: tree-twin-fit (eval = FALSE)
###################################################
## res.tree2 <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="tree", tuning=TRUE, cv=TRUE, mstop=1000, twin=TRUE, mstop2=100, 
## n.cores=2, rng=123)
## ###Number of genes selected with tree, only list 10
## sum(res.tree2$xselect==1)
## colnames(res.tree2$x)[res.tree2$xselect==1][1:10]
## pred.bj <- predict(res.tree2, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 22: dlbcl.Rnw:306-310 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 23: tree4-fit (eval = FALSE)
###################################################
## res.tree4 <- bujar(y=log(chop2[,1]), cens=chop2[,2], x=chop2[,-(1:2)], 
## learner="tree",degree=4, tuning=TRUE, cv=TRUE, mstop=100, rel.inf=TRUE, 
## n.cores=2,rng=123)
## ###Number of genes selected with tree, only list 10
## sum(res.tree4$xselect==1)
## colnames(res.tree4$x)[res.tree4$xselect==1][1:10]
## pred.bj <- predict(res.tree4, newx=rchop2[,-(1:2)])
## pred.bj <- exp(pred.bj) - 1
## group <- cut(pred.bj, breaks=c(-1, cutyear, 100), labels=c("high", "low"))
## dat.km <- data.frame(survtime=rchop$survtime, status = rchop$status, group=group)
## fit.diff <- survdiff(Surv(survtime, status) ~ group, data=dat.km)
## fit.diff


###################################################
### code chunk number 24: dlbcl.Rnw:333-337 (eval = FALSE)
###################################################
## fit.surv <- survfit(Surv(survtime, status) ~ group, data=dat.km )
## plot(fit.surv, xlab="Year past therapy",ylab="Survival probability", 
## lty = 1:2, col=c("red","blue"))
## legend(1, .1, c("High risk", "Low risk"), lty = 1:2, col=c("red","blue"))


###################################################
### code chunk number 25: dlbcl.Rnw:349-354 (eval = FALSE)
###################################################
## gene <- c("x1558999_x_at", "x212713_at", "x224043_s_at", "x229839_at", 
## "x237515_at", "x237797_at", "x242758_x_at", "x244346_at")
## par(mfrow=c(2,4))
## for(i in 1:length(gene))
## plot(res.tree4$res.fit, i.var=which(colnames(res.tree4$x) == gene[i]))  


###################################################
### code chunk number 26: dlbcl.Rnw:362-364 (eval = FALSE)
###################################################
## for(i in 1:6)
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[i,c(1, 3)]))


###################################################
### code chunk number 27: int1 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[1,c(1, 3)]))


###################################################
### code chunk number 28: int2 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[2,c(1, 3)]))


###################################################
### code chunk number 29: int3 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[3,c(1, 3)]))


###################################################
### code chunk number 30: int4 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[4,c(1, 3)]))


###################################################
### code chunk number 31: int5 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[5,c(1, 3)]))


###################################################
### code chunk number 32: int6 (eval = FALSE)
###################################################
## plot(res.tree4$res.fit, i.var=unlist(res.tree4$interactions$rank.list[6,c(1, 3)]))


