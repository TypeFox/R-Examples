### R code from vignette source 'LinRegInteractive.Rnw'

###################################################
### code chunk number 1: options
###################################################
options(prompt = " ", continue = "     ", digits = 4, show.signif.stars = FALSE)


###################################################
### code chunk number 2: LinRegInteractive.Rnw:83-86 (eval = FALSE)
###################################################
## data("creditdata")
## model.2.fac <- glm(credit ~ amount + I(amount^2)  + age + duration*teleph  
## + housing, family = binomial(link="probit"), data = creditdata)


###################################################
### code chunk number 3: LinRegInteractive.Rnw:92-93 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac) 


###################################################
### code chunk number 4: LinRegInteractive.Rnw:99-100 (eval = FALSE)
###################################################
## options(contrasts=c("contr.treatment","contr.treatment"))


###################################################
### code chunk number 5: LinRegInteractive.Rnw:106-107 (eval = FALSE)
###################################################
## options(device = "x11")


###################################################
### code chunk number 6: LinRegInteractive.Rnw:146-147 (eval = FALSE)
###################################################
## demo(VignetteFigures, package = "LinRegInteractive", ask = FALSE)


###################################################
### code chunk number 7: LinRegInteractive.Rnw:166-170 (eval = FALSE)
###################################################
## require("splines")
## model.2.fac.npamount <- glm(credit ~ bs(amount) + age + duration*teleph  
## + housing, family = binomial(link="probit"), data = creditdata)
## fxInteractive(model.2.fac.npamount) 


###################################################
### code chunk number 8: LinRegInteractive.Rnw:183-187 (eval = FALSE)
###################################################
## require("mgcv") 
## model.2.fac.mgcv <- gam(credit ~ s(amount) + age + duration*teleph + housing,  
## family = binomial(link="probit"), data = creditdata)
## fxInteractive(model.2.fac.mgcv)


###################################################
### code chunk number 9: LinRegInteractive.Rnw:220-225 (eval = FALSE)
###################################################
## require("AER")
## data("MurderRates")
## model <- glm(I(executions > 0) ~ time + income + noncauc + lfp + southern, 
## data = MurderRates, family = binomial)
## fxInteractive(model)


###################################################
### code chunk number 10: LinRegInteractive.Rnw:315-317 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## initial.values = list(amount=5000, duration=24, age=30))


###################################################
### code chunk number 11: LinRegInteractive.Rnw:326-327 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, preselect.var = "duration")


###################################################
### code chunk number 12: LinRegInteractive.Rnw:334-335 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, preselect.type = "response")


###################################################
### code chunk number 13: LinRegInteractive.Rnw:342-343 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, preselect.groups = c(1:3))


###################################################
### code chunk number 14: LinRegInteractive.Rnw:353-357 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## preselect.var     = "duration",
## snapshot.plot     = TRUE,
## graphics.filename = "D:/Temp/fig-credprobit-duration")


###################################################
### code chunk number 15: LinRegInteractive.Rnw:363-364 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, legend.width.factor = 1.1)


###################################################
### code chunk number 16: LinRegInteractive.Rnw:371-379 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## initial.values      = list(amount=5000, duration=24, age=30), 
## preselect.var       = "duration",
## preselect.type      = "marginal",
## preselect.groups    = c(2,3,5,6),
## autosave.plot       = TRUE,
## graphics.filename   = "fig-credprobit-duration-marg",
## legend.width.factor = 1.05)


###################################################
### code chunk number 17: LinRegInteractive.Rnw:388-391 (eval = FALSE)
###################################################
## data("creditdata")
## model.2.fac <- glm(credit ~ amount + I(amount^2)  + age + duration*teleph
## + housing, family = binomial(link="probit"), data = creditdata)


###################################################
### code chunk number 18: LinRegInteractive.Rnw:396-399 (eval = FALSE)
###################################################
## data("creditdata")
## model.3.fac <- glm(credit ~ amount + I(amount^2)  + age + duration*teleph 
## + housing + job, family = binomial(link="probit"), data = creditdata)


###################################################
### code chunk number 19: LinRegInteractive.Rnw:408-411 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## factor.sep = "-",
## level.sep  = ">")


###################################################
### code chunk number 20: LinRegInteractive.Rnw:419-425 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## latex2console       = TRUE,
## xtable.big.mark     = ".",
## xtable.decimal.mark = ",",
## xtable.digits       = 5,
## xtable.booktabs     = TRUE)


###################################################
### code chunk number 21: LinRegInteractive.Rnw:449-454 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## dev.height       = 11,
## dev.width        = 11,
## dev.width.legend = 5,
## dev.pointsize    = 8) 


###################################################
### code chunk number 22: LinRegInteractive.Rnw:463-467 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## preselect.var  = "amount",
## preselect.type = "response",
## ylim           = c(0,1))


###################################################
### code chunk number 23: LinRegInteractive.Rnw:475-478 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## col = c(1,"blue",2),
## lwd = rep(c(1,2),each=3))


###################################################
### code chunk number 24: LinRegInteractive.Rnw:484-489 (eval = FALSE)
###################################################
## fxInteractive(model.3.fac,  
## col = rep(c(1,"blue",2), each=4),  
## lty = c(1,2,3,4),
## lwd = rep(c(1,2), each=12),
## dev.width.legend = 8)


###################################################
### code chunk number 25: LinRegInteractive.Rnw:497-503 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## preselect.var  = "duration",
## preselect.type = "response",
## main           = "Interaction between 'duration' and factor 'teleph'",
## xlab           = "duration (months)",
## ylab           = "probability of credit default")


###################################################
### code chunk number 26: LinRegInteractive.Rnw:511-512 (eval = FALSE)
###################################################
## fxInteractive(model.3.fac, dev.width.legend = 8)


###################################################
### code chunk number 27: LinRegInteractive.Rnw:517-518 (eval = FALSE)
###################################################
## fxInteractive(model.3.fac, legend.cex = 0.7)


###################################################
### code chunk number 28: LinRegInteractive.Rnw:523-524 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, legend.pos = "top")


###################################################
### code chunk number 29: LinRegInteractive.Rnw:531-532 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, legend.add = FALSE)


###################################################
### code chunk number 30: LinRegInteractive.Rnw:537-540 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## legend.add   = FALSE,
## legend.space = TRUE)


###################################################
### code chunk number 31: LinRegInteractive.Rnw:546-547 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, legend.only = TRUE)


###################################################
### code chunk number 32: LinRegInteractive.Rnw:555-556 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, rug.ticksize = NA)


###################################################
### code chunk number 33: LinRegInteractive.Rnw:561-562 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, rug.col = "gray50")


###################################################
### code chunk number 34: LinRegInteractive.Rnw:571-572 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, vline.actual = FALSE)


###################################################
### code chunk number 35: LinRegInteractive.Rnw:578-579 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, pos.hlines = c(NA,0.56,0)) 


###################################################
### code chunk number 36: LinRegInteractive.Rnw:590-603 (eval = FALSE)
###################################################
## windows(10,7, pointsize = 10)
## layoutmatrix <- matrix(c(1,2,2), 1, 3)
## layout(layoutmatrix)
## palette(c("darkred","red","salmon","darkblue","blue","lightblue"))
## par(cex = 1, mar = c(5,5,2,2)+0.1)
## 
## fxInteractive(model.2.fac,
## preselect.var       = "amount",
## preselect.type      = "response",
## dev.defined         = TRUE,
## ylim                = c(0,1),
## legend.width.factor = 1.1,
## snapshot.plot       = TRUE)


###################################################
### code chunk number 37: LinRegInteractive.Rnw:611-616 (eval = FALSE)
###################################################
## fxInteractive(model.3.fac,
## box.type.height           = 90, 
## box.group.character.width = 6, 
## box.group.line.height     = 25, 
## dist.obj.height           = 2)


###################################################
### code chunk number 38: LinRegInteractive.Rnw:625-633 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac,
## panel.title      = "Probit Modell",
## label.button     = "Schnappschuss",
## label.slider.act = "Dargestellte Variable: ",
## label.box.type   = "Typ",
## label.types      = c("Linearer Praediktor", "Wahrscheinlichkeit",
##                      "Marginaler Effekt"),
## label.box.groups = "Gruppen")


###################################################
### code chunk number 39: LinRegInteractive.Rnw:643-649 (eval = FALSE)
###################################################
## data("munichrent03")
## require("splines")
## model.rent <- lm(rent ~ bs(yearc) + area*location + upkitchen,
## data=munichrent03)
## model.rent$data <- munichrent03
## fxInteractive(model.rent)


###################################################
### code chunk number 40: LinRegInteractive.Rnw:725-736 (eval = FALSE)
###################################################
## model.cd.manygroups <- glm(credit ~ amount + I(amount^2) + age 
## + duration*teleph + housing + intuse, family=binomial, data=creditdata)
## 
## factor.combs       <- factorCombinations(creditdata[,c("teleph",
## "housing","intuse")])
## logic.index.groups <- factor.combs$counts > 10
## index.groups       <- seq(along=factor.combs$counts)[logic.index.groups]
## 
## fxInteractive(model.cd.manygroups,
## preselect.var    = "amount",
## preselect.groups = index.groups)


###################################################
### code chunk number 41: LinRegInteractive.Rnw:741-751 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## preselect.var  = "amount",
## preselect.type = "response",
## ylim           = c(0,1),
## rug.ticksize   = 0,
## legend.width.factor = 1.1)
## segments(creditdata$amount, par("usr")[3], creditdata$amount, 
## par("fig")[3], col = rgb(0,0,0,0.2))
## # savePlot(filename = "creditdefault-customrug", type = "pdf") # Windows
## # dev.copy2pdf(file = "creditdefault-customrug.pdf") # Non-Windows


###################################################
### code chunk number 42: LinRegInteractive.Rnw:758-770 (eval = FALSE)
###################################################
## model.cd.manygroups <- glm(credit ~ amount + I(amount^2) + age
## + duration*teleph + housing + intuse, family=binomial, data=creditdata)
## index.groups <- c(1,11,21,31,41,51)
## vec.col <- NULL
## vec.col[index.groups] <- c(1:6)
## vec.lty <- NULL
## vec.lty[index.groups] <- rep(c(1,2), each = 3)
## fxInteractive(model.cd.manygroups,
## preselect.var    = "amount",
## preselect.groups = index.groups,
## col              = vec.col,
## lty              = vec.lty)


###################################################
### code chunk number 43: LinRegInteractive.Rnw:777-784 (eval = FALSE)
###################################################
## fxInteractive(model.2.fac, 
## preselect.var  = "duration",
## preselect.type = "link",
## main           = "Interaction between 'duration' and factor 'teleph'",
## main.line      = 0.5,
## cex.main       = 1,
## font.main      = 1)


