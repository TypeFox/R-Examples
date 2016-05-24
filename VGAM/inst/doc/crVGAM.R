### R code from vignette source 'crVGAM.Rnw'

###################################################
### code chunk number 1: crVGAM.Rnw:105-111
###################################################
library("VGAM")
library("VGAMdata")
ps.options(pointsize = 12)
options(width = 72, digits = 4)
options(SweaveHooks = list(fig = function() par(las = 1)))
options(prompt = "R> ", continue = "+")


###################################################
### code chunk number 2: example-posber (eval = FALSE)
###################################################
## vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
##            family = posbernoulli.t, data = pdata)


###################################################
### code chunk number 3: poz-args-posbinomial
###################################################
args(posbinomial)


###################################################
### code chunk number 4: poz-args-posbernoulli-t
###################################################
args(posbernoulli.t)


###################################################
### code chunk number 5: poz-args-posbernoulli-b
###################################################
args(posbernoulli.b)


###################################################
### code chunk number 6: poz-args-posbernoulli-tb
###################################################
args(posbernoulli.tb)


###################################################
### code chunk number 7: poz-posbernoulli-tb-gen (eval = FALSE)
###################################################
## vglm(..., family = posbernoulli.tb(parallel.b = TRUE ~ 0, parallel.t = TRUE ~ 0,
##                                    drop.b = TRUE ~ 0))


###################################################
### code chunk number 8: eg-deermice-look
###################################################
head(deermice, 4)


###################################################
### code chunk number 9: example1-model
###################################################
deermice <- within(deermice, {
  age <- 2 - as.numeric(age) 
  sex <- 1 - as.numeric(sex)
})


###################################################
### code chunk number 10: example2-model
###################################################
M.0 <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1,
            posbernoulli.t(parallel = TRUE ~ 1), data = deermice)
M.b <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
            posbernoulli.b, data = deermice)
M.t <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
            posbernoulli.t, data = deermice)
M.h <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
            posbernoulli.t(parallel = TRUE ~ weight + sex + age), data = deermice)
M.th <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
             posbernoulli.t, data = deermice)
M.tb <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ 1, 
             posbernoulli.tb, data = deermice)
M.bh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
             posbernoulli.b, data = deermice)
M.tbh <- vglm(cbind(y1, y2, y3, y4, y5, y6) ~ weight + sex + age,
              posbernoulli.tb, data = deermice)


###################################################
### code chunk number 11: eg-deermice-Nhat
###################################################
c(M.bh@extra$N.hat, M.bh@extra$SE.N.hat)
c(logLik(M.bh), AIC(M.bh))


###################################################
### code chunk number 12: maketable
###################################################

Table <- rbind(c(round(M.tbh@extra$N.hat,2), 
                 round(M.bh@extra$N.hat,2), 
                 round(M.tb@extra$N.hat,2), 
                 round(M.th@extra$N.hat,2), 
                 round(M.h@extra$N.hat,2), 
                 round(M.b@extra$N.hat,2),
                 round(M.t@extra$N.hat,2), 
                 round(M.0@extra$N.hat,2)), 
               
              c(round(M.tbh@extra$SE.N.hat,2), 
                round(M.bh@extra$SE.N.hat,2), 
                round(M.tb@extra$SE.N.hat,2),
                round(M.th@extra$SE.N.hat,2), 
                round(M.h@extra$SE.N.hat,2), 
                round(M.b@extra$SE.N.hat,2),
                round(M.t@extra$SE.N.hat,2), 
                round(M.0@extra$SE.N.hat,2)), 
               
              -2*c(round(logLik(M.tbh),2), 
                   round(logLik(M.bh),2), 
                   round(logLik(M.tb),2), 
                   round(logLik(M.th),2), 
                   round(logLik(M.h),2), 
                   round(logLik(M.b),2),
                   round(logLik(M.t),2), 
                   round(logLik(M.0),2)), 
               
              c(round(AIC(M.tbh),2), 
                round(AIC(M.bh),2), 
                round(AIC(M.tb),2), 
                round(AIC(M.th),2),
                round(AIC(M.h),2), 
                round(AIC(M.b),2), 
                round(AIC(M.t),2), 
                round(AIC(M.0),2)));

colnames(Table) <- c("M.tbh", "M.bh", "M.tb", 
                     "M.th", "M.h", "M.b", "M.t", "M.0");
rownames(Table) <- c("N.hat", "SE","-2ln(L)", "AIC");


###################################################
### code chunk number 13: example2-table
###################################################
Table


###################################################
### code chunk number 14: poz-posbernoulli-eg-deermice-coefs
###################################################
round(coef(M.bh), 2)
round(sqrt(diag(vcov(M.bh))), 2)


###################################################
### code chunk number 15: poz-posbernoulli-eg-deermice-smooth
###################################################
fit.bh <- vgam(cbind(y1, y2, y3, y4, y5, y6) ~ s(weight, df = 3) + sex + age,
               posbernoulli.b, data = deermice)
plot(fit.bh, se = TRUE, las = 1, lcol = "blue", scol = "orange",
     rcol = "purple", scale = 5)


###################################################
### code chunk number 16: poz-posbernoulli-eg-deermice-summary
###################################################
summary(fit.bh)


###################################################
### code chunk number 17: poz-posbernoulli-eg-deermice-smooth-shadow (eval = FALSE)
###################################################
## plot(fit.bh, se = TRUE, las = 1, lcol = "blue", scol = "orange",
##      rcol = "purple", scale = 5, mgp = c(2.0, 1, 0))


###################################################
### code chunk number 18: plot-deermice
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow = c(2, 2))
par(las = 1, cex = 1.1, mar = c(3.8, 4, 0.5, 0.2) + 0.1)
par(mgp = c(2.3, 1, 0))  # Default is c(3, 1, 0)



plot(fit.bh, se = TRUE, las = 1, lcol = "blue", scol = "orange",
     rcol = "purple", scale = 5, mgp = c(2.0, 1, 0))

# < < poz-posbernoulli-eg-deermice-smooth-shadow> >




###################################################
### code chunk number 19: birds91read
###################################################
data("prinia", package = "VGAM")


###################################################
### code chunk number 20: example2a
###################################################
head(prinia, 4)[, 1:4]


###################################################
### code chunk number 21: example2b
###################################################
M.h.GAM <- 
  vgam(cbind(cap, noncap) ~ s(length, df = 3) + fat, 
       posbinomial(omit.constant = TRUE, parallel = TRUE ~ s(length, df = 3) + fat),
       data = prinia)
M.h.GAM@extra$N.hat     
M.h.GAM@extra$SE.N.hat  


###################################################
### code chunk number 22: eg-bird-smooth-shadow1
###################################################
plot.info <- plot(M.h.GAM,
                  se = TRUE, las = 1, plot.arg = FALSE,
                  lcol = "blue",
                  scol = "orange",
                  rcol = "purple",
                  scale = 5)


###################################################
### code chunk number 23: eg-bird-smooth-shadow2 (eval = FALSE)
###################################################
## info.fit2 <- plot.info@preplot[[1]]
## fat.effect <- coef(M.h.GAM)["fat"] 
## intercept <- coef(M.h.GAM)["(Intercept)"]  
## 
## ooo <- order(info.fit2$x)
## centering.const <- mean(prinia$length) - coef(M.h.GAM)["s(length, df = 3)"]
## 
## plotframe <- data.frame(lin.pred.b = intercept + fat.effect * 1 +
##                                      centering.const + info.fit2$y[ooo],
##                         lin.pred.0 = intercept + fat.effect * 0 +
##                                      centering.const + info.fit2$y[ooo],
##                         x2 = info.fit2$x[ooo])
## 
## plotframe <- transform(plotframe,
##                        up.lin.pred.b = lin.pred.b + 2*info.fit2$se.y[ooo],
##                        lo.lin.pred.b = lin.pred.b - 2*info.fit2$se.y[ooo],
##                        up.lin.pred.0 = lin.pred.0 + 2*info.fit2$se.y[ooo],
##                        lo.lin.pred.0 = lin.pred.0 - 2*info.fit2$se.y[ooo])
## 
## plotframe <- transform(plotframe,
##                        fv.b    = logit(lin.pred.b,    inverse = TRUE),
##                        up.fv.b = logit(up.lin.pred.b, inverse = TRUE),
##                        lo.fv.b = logit(lo.lin.pred.b, inverse = TRUE),
##                        fv.0    = logit(lin.pred.0,    inverse = TRUE),
##                        up.fv.0 = logit(up.lin.pred.0, inverse = TRUE),
##                        lo.fv.0 = logit(lo.lin.pred.0, inverse = TRUE))
## 
## with(plotframe,
##      matplot(x2, cbind(up.fv.b, fv.b, lo.fv.b), type = "l", col = "blue",
##              lty = c(2, 1, 2), las = 1, cex.lab = 1.5, lwd = 2,
##              main = "", ylab = "", xlab = "Wing length (standardized)"))
## mtext( ~ hat(p), side = 2, cex = 1.4, line = 4, adj = 0.5, las = 1)
## with(plotframe, matlines(x2, cbind(up.fv.0, fv.0, lo.fv.0),
##                          col = "darkorange", lty = c(2, 1, 2)), lwd = 2)
## legend("topleft", legend = c("Fat present", "Fat not present"), bty = "n",
##        lwd = 2, col = c("blue", "darkorange"), merge = TRUE, cex = 1.5)


###################################################
### code chunk number 24: plot-bird
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfrow = c(1, 1))



info.fit2 <- plot.info@preplot[[1]]
fat.effect <- coef(M.h.GAM)["fat"] 
intercept <- coef(M.h.GAM)["(Intercept)"]  

ooo <- order(info.fit2$x)
centering.const <- mean(prinia$length) - coef(M.h.GAM)["s(length, df = 3)"]

plotframe <- data.frame(lin.pred.b = intercept + fat.effect * 1 +
                                     centering.const + info.fit2$y[ooo],
                        lin.pred.0 = intercept + fat.effect * 0 +
                                     centering.const + info.fit2$y[ooo],
                        x2 = info.fit2$x[ooo])

plotframe <- transform(plotframe,
                       up.lin.pred.b = lin.pred.b + 2*info.fit2$se.y[ooo],
                       lo.lin.pred.b = lin.pred.b - 2*info.fit2$se.y[ooo],
                       up.lin.pred.0 = lin.pred.0 + 2*info.fit2$se.y[ooo],
                       lo.lin.pred.0 = lin.pred.0 - 2*info.fit2$se.y[ooo])

plotframe <- transform(plotframe,
                       fv.b    = logit(lin.pred.b,    inverse = TRUE),
                       up.fv.b = logit(up.lin.pred.b, inverse = TRUE),
                       lo.fv.b = logit(lo.lin.pred.b, inverse = TRUE),
                       fv.0    = logit(lin.pred.0,    inverse = TRUE),
                       up.fv.0 = logit(up.lin.pred.0, inverse = TRUE),
                       lo.fv.0 = logit(lo.lin.pred.0, inverse = TRUE))

with(plotframe,
     matplot(x2, cbind(up.fv.b, fv.b, lo.fv.b), type = "l", col = "blue",
             lty = c(2, 1, 2), las = 1, cex.lab = 1.5, lwd = 2,
             main = "", ylab = "", xlab = "Wing length (standardized)"))
mtext( ~ hat(p), side = 2, cex = 1.4, line = 4, adj = 0.5, las = 1)
with(plotframe, matlines(x2, cbind(up.fv.0, fv.0, lo.fv.0),
                         col = "darkorange", lty = c(2, 1, 2)), lwd = 2)
legend("topleft", legend = c("Fat present", "Fat not present"), bty = "n",
       lwd = 2, col = c("blue", "darkorange"), merge = TRUE, cex = 1.5)



# < < eg-bird-smooth-shadow2 > >
    
    
    


###################################################
### code chunk number 25: poz-posbernoulli-tb-huggins89t1-data
###################################################
head(Huggins89table1, 4)


###################################################
### code chunk number 26: poz-posbernoulli-tb-huggins89t1-look
###################################################
Hdata <- transform(Huggins89table1, x3.tij = t01,
                   T02 = t02, T03 = t03, T04 = t04, T05 = t05, T06 = t06,
                   T07 = t07, T08 = t08, T09 = t09, T10 = t10)
Hdata <- subset(Hdata,
                y01 + y02 + y03 + y04 + y05 + y06 + y07 + y08 + y09 + y10 > 0)


###################################################
### code chunk number 27: poz-posbernoulli-th-huggins89t0-fit
###################################################
fit.th <-
   vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10) ~  x2 + x3.tij,
        xij = list(x3.tij ~ t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                            t09 + t10 - 1),
        posbernoulli.t(parallel.t = TRUE ~ x2 + x3.tij), 
        data = Hdata, trace = FALSE, 
        form2 = ~ x2 + x3.tij + t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 +
                                t09 + t10)


###################################################
### code chunk number 28: poz-posbernoulli-th-huggins89t0-constraints
###################################################
constraints(fit.th, matrix = TRUE)


###################################################
### code chunk number 29: poz-posbernoulli-tbh-huggins89t1-fit
###################################################
fit.tbh <-
  vglm(cbind(y01, y02, y03, y04, y05, y06, y07, y08, y09, y10) ~  x2 + x3.tij,
       xij = list(x3.tij ~ t01 + t02 + t03 + t04 + t05 + t06 +
                                 t07 + t08 + t09 + t10 +
                                 T02 + T03 + T04 + T05 + T06 +
                                 T07 + T08 + T09 + T10 - 1),
       posbernoulli.tb(parallel.t = TRUE ~ x2 + x3.tij),
       data = Hdata, trace = FALSE,
       form2 = ~  x2 + x3.tij +
                  t01 + t02 + t03 + t04 + t05 + t06 + t07 + t08 + t09 + t10 +
                        T02 + T03 + T04 + T05 + T06 + T07 + T08 + T09 + T10)


###################################################
### code chunk number 30: poz-posbernoulli-tbh-huggins89t1-aic
###################################################
c(logLik(fit.th), AIC(fit.th))
c(logLik(fit.tbh), AIC(fit.tbh))


###################################################
### code chunk number 31: poz-posbernoulli-tb-huggins89t1-constraints
###################################################
head(constraints(fit.tbh, matrix = TRUE), 4)
tail(constraints(fit.tbh, matrix = TRUE), 4)


###################################################
### code chunk number 32: poz-posbernoulli-tb-huggins89t1-coefs
###################################################
coef(fit.tbh)
sqrt(diag(vcov(fit.tbh))) 


###################################################
### code chunk number 33: poz-posbernoulli-tb-huggins89t1-Nhat
###################################################
fit.tbh@extra$N.hat   
fit.tbh@extra$SE.N.hat


###################################################
### code chunk number 34: poz-posbernoulli-tbh-huggins89t1-fit-Select
###################################################
Hdata <- subset(Huggins89table1, rowSums(Select(Huggins89table1, "y")) > 0)
Hdata.T <- Select(Hdata, "t")
colnames(Hdata.T) <- gsub("t", "T", colnames(Hdata.T))
Hdata <- data.frame(Hdata, Hdata.T)
Hdata <- transform(Hdata, x3.tij = y01)
Form2 <- Select(Hdata, prefix = TRUE, as.formula = TRUE)
Xij   <- Select(Hdata, c("t", "T"), as.formula = TRUE,
                sort = FALSE, rhs = "0", lhs = "x3.tij", exclude = "T01")
fit.tbh <- vglm(Select(Hdata, "y") ~ x2 + x3.tij,
                form2 = Form2,  xij = list(Xij),
                posbernoulli.tb(parallel.t = TRUE ~ x2 + x3.tij),
                data = Hdata, trace = FALSE)
coef(fit.tbh)


###################################################
### code chunk number 35: poz-posbernoulli-bh-ephemeral-method1
###################################################
deermice <- transform(deermice, Lag1 = y1)
M.tbh.lag1 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag1,
       posbernoulli.tb(parallel.t = FALSE ~ 0,
                       parallel.b = FALSE ~ 0,
                       drop.b = FALSE ~ 1),
       xij = list(Lag1 ~ fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                         fill(y5) + fill(y6) +
                         y1 + y2 + y3 + y4 + y5),
       form2 = ~ sex + weight + Lag1 +
                 fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                 fill(y5) + fill(y6) +
                 y1 + y2 + y3 + y4 + y5 + y6,
       data = deermice)
coef(M.tbh.lag1)


###################################################
### code chunk number 36: poz-posbernoulli-bh-ephemeral-method2
###################################################
deermice <- transform(deermice, Lag1 = y1)
deermice <- transform(deermice, f1 = y1, f2 = y1, f3 = y1, f4 = y1,
                                f5 = y1, f6 = y1)
tau <- 6
H2 <- H3 <- cbind(rep(1, 2*tau-1))
H4 <- cbind(c(rep(0, tau), rep(1, tau-1)))
M.tbh.lag1.method2 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag1,
       posbernoulli.tb(parallel.b = TRUE ~ 0, parallel.t = TRUE ~ 0),
       constraints = list("(Intercept)" = cbind(H4, 1), sex = H2, weight= H3, 
                          Lag1 = H4),
       xij = list(Lag1 ~ f1 + f2 + f3 + f4 + f5 + f6 +
                         y1 + y2 + y3 + y4 + y5),
       form2 = Select(deermice, prefix = TRUE, as.formula = TRUE),
       data = deermice)
coef(M.tbh.lag1.method2)


###################################################
### code chunk number 37: poz-posbernoulli-bh-ephemeral-lag2
###################################################
deermice <- transform(deermice, Lag2 = y1)
M.bh.lag2 <-
  vglm(cbind(y1, y2, y3, y4, y5, y6) ~ sex + weight + Lag2,
       posbernoulli.tb(parallel.t = FALSE ~ 0,
                       parallel.b = FALSE ~ 0,
                       drop.b = FALSE ~ 1),
       xij = list(Lag2 ~ fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                         fill(y5) + fill(y6) +
                         y1 + pmax(y1, y2) + pmax(y2, y3) + pmax(y3, y4) + 
                         pmax(y4, y5)),
       form2 = ~ sex + weight + Lag2 +
                 fill(y1) + fill(y2) + fill(y3) + fill(y4) +
                 fill(y5) + fill(y6) +
                 y1 + pmax(y1, y2) + pmax(y2, y3) + pmax(y3, y4) + 
                 pmax(y4, y5) + y6,
       data = deermice)
coef(M.bh.lag2)


