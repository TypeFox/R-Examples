### R code from vignette source 'loglinear-birth.Rnw'

###################################################
### code chunk number 1: loglinear-birth.Rnw:13-14 (eval = FALSE)
###################################################
## options(width=80)


###################################################
### code chunk number 2: loglinear-birth.Rnw:17-21 (eval = FALSE)
###################################################
## library(catdata)
## 
## data(birth)
## attach(birth)


###################################################
### code chunk number 3: loglinear-birth.Rnw:27-30 (eval = FALSE)
###################################################
## table1 <- table(Sex, Membranes, Cesarean, Induced)
## 
## ftable(table1)


###################################################
### code chunk number 4: loglinear-birth.Rnw:38-50 (eval = FALSE)
###################################################
## m4 <- loglin(table1, margin=list(c(1,2,3,4)), fit=TRUE)
## cat("deviance(m4)=", m4$lrt, "df(m4)=", m4$df, "\n")
## 
## m3 <- loglin(table1, margin=list(c(1,2,3), c(1,2,4), c(1,3,4), c(2,3,4)), fit=TRUE)
## cat("deviance(m3)=", m3$lrt, "df(m3)=", m3$df, "\n")
## 
## m2 <- loglin(table1, margin=list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)), 
##              fit=TRUE)
## cat("deviance(m2)=", m2$lrt, "df(m2)=", m2$df, "\n")
## 
## m1 <- loglin(table1, margin=list(c(1), c(2), c(3), c(4)), fit=TRUE)
## cat("deviance(m1)=", m1$lrt, "df(m1)=", m1$df, "\n")


###################################################
### code chunk number 5: loglinear-birth.Rnw:55-67 (eval = FALSE)
###################################################
## (df34 <- m3$df - m4$df)
## (dev34 <- m3$lrt - m4$lrt)
## 1-pchisq(dev34, df34)
## 
## 
## (df23 <- m2$df - m3$df)
## (dev23 <- m2$lrt - m3$lrt)
## 1-pchisq(dev23, df23)
## 
## (df12 <- m1$df - m2$df)
## (dev12 <- m1$lrt - m2$lrt)
## 1-pchisq(dev12, df12)


###################################################
### code chunk number 6: loglinear-birth.Rnw:72-95 (eval = FALSE)
###################################################
## m2.GM <- loglin(table1, margin=list(c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)), 
##                 fit=TRUE)
## cat("deviance(m2.GM)=", m2.GM$lrt, "df(m2.GM)=", m2.GM$df, "\n")
## 
## m2.MC <- loglin(table1, margin=list(c(1,2), c(1,3), c(1,4), c(2,4), c(3,4)), 
##                 fit=TRUE)
## cat("deviance(m2.MC)=", m2.MC$lrt, "df(m2.MC)=", m2.MC$df, "\n")
## 
## m2.CI <- loglin(table1, margin=list(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4)), 
##                 fit=TRUE)
## cat("deviance(m2.CI)=", m2.CI$lrt, "df(m2.CI)=", m2.CI$df, "\n")
## 
## m2.GI <- loglin(table1, margin=list(c(1,2), c(1,3), c(2,3), c(2,4), c(3,4)), 
##                 fit=TRUE)
## cat("deviance(m2.GI)=", m2.GI$lrt, "df(m2.GI)=", m2.GI$df, "\n")
## 
## m2.GC <- loglin(table1, margin=list(c(1,2), c(1,4), c(2,3), c(2,4), c(3,4)), 
##                 fit=TRUE)
## cat("deviance(m2.GC)=", m2.GC$lrt, "df(m2.GC)=", m2.GC$df, "\n")
## 
## m2.MI <- loglin(table1, margin=list(c(1,2), c(1,3), c(1,4), c(2,3), c(3,4)), 
##                 fit=TRUE)
## cat("deviance(m2.MI)=", m2.MI$lrt, "df(m2.MI)=", m2.MI$df, "\n")


###################################################
### code chunk number 7: loglinear-birth.Rnw:101-112 (eval = FALSE)
###################################################
## 1 - pchisq(m2.GM$lrt - m2$lrt, 1)
## 
## 1 - pchisq(m2.MC$lrt - m2$lrt, 1)
## 
## 1 - pchisq(m2.CI$lrt - m2$lrt, 1)
## 
## 1 - pchisq(m2.GI$lrt - m2$lrt, 1)
## 
## 1 - pchisq(m2.GC$lrt - m2$lrt, 1)
## 
## 1 - pchisq(m2.MI$lrt - m2$lrt, 1)


###################################################
### code chunk number 8: loglinear-birth.Rnw:118-123 (eval = FALSE)
###################################################
## m2.GM.GI.GC<- loglin(table1, margin=list(c(1), c(2,3), c(2,4), c(3,4)), fit=TRUE)
## cat("deviance(m2.GM.GI.GC)=", m2.GM.GI.GC$lrt, "df(m2.GM.GI.GC)=", m2.GM.GI.GC$df, 
##     "\n")
## 
## 1 - pchisq(m2.GM.GI.GC$lrt - m2$lrt, m2.GM.GI.GC$df - m2$df)


###################################################
### code chunk number 9: loglinear-birth.Rnw:129-133 (eval = FALSE)
###################################################
## m2.G<- loglin(table1, margin=list(c(2,3), c(2,4), c(3,4)), fit=TRUE)
## cat("deviance(m2.G)=", m2.G$lrt, "df(m2.G)=", m2.G$df, "\n")
## 
## 1 - pchisq(m2.G$lrt - m2$lrt, m2.G$df - m2$df)


