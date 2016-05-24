### R code from vignette source 'twang.rnw'

###################################################
### code chunk number 1: twang.rnw:104-105
###################################################
options(width=60)


###################################################
### code chunk number 2: twang.rnw:124-126
###################################################
library(twang)
set.seed(1)


###################################################
### code chunk number 3: twang.rnw:134-135
###################################################
data(lalonde)


###################################################
### code chunk number 4: twang.rnw:158-168
###################################################
ps.lalonde <- ps(treat ~ age + educ + black + hispan + nodegree +
                         married + re74 + re75,
                 data = lalonde,
                 n.trees=5000,
                 interaction.depth=2,
                 shrinkage=0.01,
                 perm.test.iters=0,
                 stop.method=c("es.mean","ks.max"),                 
                 estimand = "ATT",
                 verbose=FALSE)


###################################################
### code chunk number 5: iterPt
###################################################
    plot(ps.lalonde)


###################################################
### code chunk number 6: iterPt2
###################################################
    plot(ps.lalonde, subset = 2)


###################################################
### code chunk number 7: twang.rnw:366-369
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees,
        plot=FALSE)


###################################################
### code chunk number 8: twang.rnw:375-377
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees)


###################################################
### code chunk number 9: twang.rnw:406-407
###################################################
options(width=85)


###################################################
### code chunk number 10: twang.rnw:410-412
###################################################
lalonde.balance <- bal.table(ps.lalonde)
lalonde.balance


###################################################
### code chunk number 11: twang.rnw:415-416
###################################################
options(width=60)


###################################################
### code chunk number 12: twang.rnw:489-498
###################################################
library(xtable)
pretty.tab <- lalonde.balance$ks.max.ATT[,c("tx.mn","ct.mn","ks")]
pretty.tab <- cbind(pretty.tab, lalonde.balance$unw[,"ct.mn"])
names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","E(Y0|t=0)")
xtable(pretty.tab,
       caption = "Balance of the treatment and comparison groups",
       label = "tab:balance",
       digits = c(0, 2, 2, 2, 2),
       align=c("l","r","r","r","r"))


###################################################
### code chunk number 13: twang.rnw:512-513
###################################################
summary(ps.lalonde)


###################################################
### code chunk number 14: twang.rnw:589-590
###################################################
plot(ps.lalonde, plots=2)


###################################################
### code chunk number 15: twang.rnw:635-636
###################################################
plot(ps.lalonde, plots=3)


###################################################
### code chunk number 16: twang.rnw:651-652
###################################################
plot(ps.lalonde, plots = 4)


###################################################
### code chunk number 17: twang.rnw:667-668
###################################################
plot(ps.lalonde, plots = 5)


###################################################
### code chunk number 18: twang.rnw:703-704
###################################################
plot(ps.lalonde, plots = 3, subset = 2)


###################################################
### code chunk number 19: twang.rnw:718-719
###################################################
library(survey)


###################################################
### code chunk number 20: twang.rnw:731-733
###################################################
lalonde$w <- get.weights(ps.lalonde, stop.method="es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=lalonde)


###################################################
### code chunk number 21: twang.rnw:736-736
###################################################



###################################################
### code chunk number 22: twang.rnw:770-772
###################################################
glm1 <- svyglm(re78 ~ treat, design=design.ps)
summary(glm1)


###################################################
### code chunk number 23: twang.rnw:804-806
###################################################
glm2 <- svyglm(re78 ~ treat + nodegree, design=design.ps)
summary(glm2)


###################################################
### code chunk number 24: twang.rnw:818-822
###################################################
glm3 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.ps)
summary(glm3)


###################################################
### code chunk number 25: twang.rnw:831-835
###################################################
glm4 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree +
                  married + re74 + re75,
           data=lalonde)
summary(glm4)


###################################################
### code chunk number 26: twang.rnw:838-841
###################################################
glm5 <- lm(sqrt(re78) ~ treat + age + educ + black + hispan + nodegree +
                        married + sqrt(re74) + sqrt(re75),
           data=lalonde)


###################################################
### code chunk number 27: twang.rnw:867-873
###################################################
ps.logit <- glm(treat ~ age + educ + black + hispan + nodegree +
                        married + re74 + re75,
                data = lalonde,
                family = binomial)
lalonde$w.logit <- rep(1,nrow(lalonde))
lalonde$w.logit[lalonde$treat==0] <- exp(predict(ps.logit,subset(lalonde,treat==0)))


###################################################
### code chunk number 28: twang.rnw:887-894
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 29: twang.rnw:902-909
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 30: twang.rnw:921-929
###################################################
pretty.tab <- bal.table(bal.logit)[[2]][,c("tx.mn","ct.mn","ks")]
pretty.tab <- cbind(pretty.tab, bal.table(bal.logit)[[1]]$ct.mn)
names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","E(Y0|t=0)")
xtable(pretty.tab,
       caption = "Logistic regression estimates of the propensity scores",
       label = "tab:balancelogit",
       digits = c(0, 2, 2, 2, 2),
       align=c("l","r","r","r","r"))


###################################################
### code chunk number 31: twang.rnw:939-956
###################################################
bal.gbm <- dx.wts(ps.lalonde,
                  data=lalonde, estimand = "ATE",
                  vars=c("age","educ","black","hispan","nodegree","married","re74","re75"),
                  treat.var="treat",
                  perm.test.iters=0)
pretty.tab <- rbind(bal.logit$summary.tab,
                    bal.gbm$summary.tab[-1,])
rownames(pretty.tab) <- pretty.tab$type
rownames(pretty.tab)[2] <- "logit"
pretty.tab <- pretty.tab[,c("n.treat","ess.ctrl","max.es","mean.es","max.ks","mean.ks")]
xtable(pretty.tab,
       caption = "Summary of the balancing properties of logistic regression and gbm",
       label = "tab:balancecompare",
       digits = c(0, 0, 2, 2, 2, 2, 2),
       #digits = c(0,0, rep(2,9)),
       align=c("l","r","r","r","r","r","r"))
       #align = rep("l", 11))


###################################################
### code chunk number 32: twang.rnw:959-962
###################################################
design.logit <- svydesign(ids=~1, weights=~w.logit, data=lalonde)
glm6 <- svyglm(re78 ~ treat, design=design.logit)
summary(glm6)


###################################################
### code chunk number 33: twang.rnw:975-978
###################################################
glm7 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.logit)


###################################################
### code chunk number 34: twang.rnw:1057-1060
###################################################
data(lindner)
table(lindner$sixMonthSurvive, lindner$abcix)
chisq.test(table(lindner$sixMonthSurvive, lindner$abcix))


###################################################
### code chunk number 35: twang.rnw:1071-1075
###################################################
set.seed(1)
ps.lindner <- ps(abcix ~ stent + height + female + diabetic + 
                 acutemi + ejecfrac + ves1proc, data = lindner,
                 verbose = FALSE, estimand = "ATE")


###################################################
### code chunk number 36: twang.rnw:1087-1088
###################################################
options(width=85)


###################################################
### code chunk number 37: twang.rnw:1091-1092
###################################################
bal.table(ps.lindner)


###################################################
### code chunk number 38: twang.rnw:1095-1096
###################################################
options(width = 60)


###################################################
### code chunk number 39: twang.rnw:1118-1119
###################################################
plot(ps.lindner, plots = 1)


###################################################
### code chunk number 40: twang.rnw:1126-1127
###################################################
plot(ps.lindner, plots = 2)


###################################################
### code chunk number 41: twang.rnw:1133-1134
###################################################
plot(ps.lindner, plots = 3)


###################################################
### code chunk number 42: twang.rnw:1141-1142
###################################################
plot(ps.lindner, plots = 4)


###################################################
### code chunk number 43: twang.rnw:1148-1149
###################################################
plot(ps.lindner, plots = 5)


###################################################
### code chunk number 44: twang.rnw:1162-1163
###################################################
summary(ps.lindner)


###################################################
### code chunk number 45: twang.rnw:1170-1173
###################################################
lindner$w <- get.weights(ps.lindner, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights = ~w, data = lindner)
svychisq(~sixMonthSurvive + abcix, design = design.ps)


###################################################
### code chunk number 46: twang.rnw:1233-1234
###################################################
data(egsingle)


###################################################
### code chunk number 47: twang.rnw:1242-1243
###################################################
tmp <- tapply(egsingle$grade, egsingle$childid, unique) 


###################################################
### code chunk number 48: twang.rnw:1251-1252
###################################################
tmp <- lapply(tmp, function(x){return(x %in% 1:4)})


###################################################
### code chunk number 49: twang.rnw:1259-1260
###################################################
tmp <- lapply(tmp, sum)


###################################################
### code chunk number 50: twang.rnw:1265-1266
###################################################
tmp <- sapply(tmp, function(x){as.numeric(x == 4)})


###################################################
### code chunk number 51: twang.rnw:1271-1274
###################################################
tmp <- data.frame(tmp)
names(tmp) <- "resp"
tmp$childid <- row.names(tmp)


###################################################
### code chunk number 52: twang.rnw:1279-1280
###################################################
egsingle <- merge(egsingle, tmp)


###################################################
### code chunk number 53: twang.rnw:1287-1288
###################################################
egsingle.one <-unique(egsingle[,-c(3:6)])


###################################################
### code chunk number 54: twang.rnw:1293-1295
###################################################
egsingle.one$race <- as.factor(race <- ifelse(egsingle.one$black==1, 1,
                                         ifelse(egsingle.one$hispanic==1, 2, 3)))


###################################################
### code chunk number 55: twang.rnw:1307-1313
###################################################
egsingle.ps <-  ps(resp ~ race + female + size + lowinc + mobility,
      data=egsingle.one,
      stop.method=c("es.mean","ks.max"),
      n.trees=5000,
      verbose=FALSE,
      estimand = "ATE")


###################################################
### code chunk number 56: twang.rnw:1324-1325
###################################################
plot(egsingle.ps)


###################################################
### code chunk number 57: twang.rnw:1363-1364
###################################################
egsingle.one$wgt <- get.weights(egsingle.ps, stop.method="ks.max")


###################################################
### code chunk number 58: twang.rnw:1373-1375
###################################################
egtmp <- rbind(data.frame(egsingle.one, nr2=1, wgt2=1), 
                data.frame(egsingle.one, nr2=0, wgt2=egsingle.one$wgt)[egsingle.one$resp==1,])


###################################################
### code chunk number 59: twang.rnw:1381-1386
###################################################
egdxwts <- dx.wts(x=egtmp$wgt2, 
                   data=egtmp, 
                   estimand="ATT", 
                   vars=c("race", "female", "size",  "lowinc",  "mobility"),
                   treat.var="nr2")


###################################################
### code chunk number 60: twang.rnw:1389-1396
###################################################
pretty.tab<-bal.table(egdxwts)[[2]][,c("tx.mn","ct.mn","std.eff.sz","ks")]
names(pretty.tab) <- c("OverallS Sample","Weighted responders","Std ES","KS")
xtable(pretty.tab,
        caption = "Balance of the nonrespondents and respondents",
        label = "tab:balance2",
        digits = c(0, 2, 2, 2, 2),
        align=c("l","r","r","r","r"))


###################################################
### code chunk number 61: twang.rnw:1405-1408
###################################################
egsinge.resp <- merge(subset(egsingle, subset=resp==1),
                        subset(egsingle.one, subset=resp==1,
                               select=c(childid, wgt)) )


