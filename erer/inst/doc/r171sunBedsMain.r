# Title: R Program for Sun (2011 FPE)
library(apt); library(vars); setwd('C:/aErer')
options(width = 100, stringsAsFactors = FALSE)

# -------------------------------------------------------------------------
# 1. Data and summary statistics
# Price data for China and Vietnam are saved as 'daVich'
data(daVich); head(daVich); tail(daVich); str(daVich)
prVi <- daVich[, 1]; prCh <- daVich[, 2]
(dog <- t(bsStat(y = daVich, digits = c(3, 3)))) 
dog2 <- data.frame(item = rownames(dog), CH.level = dog[, 2], 
  CH.diff = '__', VI.level = dog[, 1], VI.diff = '__')[2:6, ]
rownames(dog2) <- 1:nrow(dog2); str(dog2); dog2

# -------------------------------------------------------------------------
# 2. Unit root test (Table 1)
ch.t1 <- ur.df(type = 'trend', lags = 3,  y = prCh); slotNames(ch.t1)
ch.d1 <- ur.df(type = 'drift', lags = 3,  y = prCh)
ch.t2 <- ur.df(type = 'trend', lags = 3,  y = diff(prCh))
ch.d2 <- ur.df(type = 'drift', lags = 3,  y = diff(prCh))
vi.t1 <- ur.df(type = 'trend', lags = 12, y = prVi)
vi.d1 <- ur.df(type = 'drift', lags = 11, y = prVi)
vi.t2 <- ur.df(type = 'trend', lags = 10, y = diff(prVi))
vi.d2 <- ur.df(type = 'drift', lags = 10, y = diff(prVi))
dog2[6, ] <- c('ADF with trend', 
  paste(round(ch.t1@teststat[1], digits = 3), '[', 3,  ']', sep = ''),
  paste(round(ch.t2@teststat[1], digits = 3), '[', 3,  ']', sep = ''),
  paste(round(vi.t1@teststat[1], digits = 3), '[', 12, ']', sep = ''),
  paste(round(vi.t2@teststat[1], digits = 3), '[', 10, ']', sep = ''))
dog2[7, ] <- c('ADF with drift', 
  paste(round(ch.d1@teststat[1], digits = 3), '[', 3,  ']', sep = ''),
  paste(round(ch.d2@teststat[1], digits = 3), '[', 3,  ']', sep = ''),
  paste(round(vi.d1@teststat[1], digits = 3), '[', 11, ']', sep = ''),
  paste(round(vi.d2@teststat[1], digits = 3), '[', 10, ']', sep = ''))
(table.1 <- dog2) 

# -------------------------------------------------------------------------
# 3. Johansen-Juselius and Engle-Granger cointegration analyses
# JJ cointegration
VARselect(daVich, lag.max = 12, type = 'const')
summary(VAR(daVich, type = 'const', p = 1))
K <- 5; two <- cbind(prVi, prCh)
summary(j1 <- ca.jo(x = two, type = 'eigen', ecdet = 'trend', K = K))
summary(j2 <- ca.jo(x = two, type = 'eigen', ecdet = 'const', K = K))
summary(j3 <- ca.jo(x = two, type = 'eigen', ecdet = 'none' , K = K))
summary(j4 <- ca.jo(x = two, type = 'trace', ecdet = 'trend', K = K))
summary(j5 <- ca.jo(x = two, type = 'trace', ecdet = 'const', K = K))
summary(j6 <- ca.jo(x = two, type = 'trace', ecdet = 'none' , K = K))
slotNames(j1)
out1 <- cbind('eigen', 'trend', K, round(j1@teststat, digits = 3), j1@cval)
out2 <- cbind('eigen', 'const', K, round(j2@teststat, digits = 3), j2@cval)
out3 <- cbind('eigen', 'none',  K, round(j3@teststat, digits = 3), j3@cval)
out4 <- cbind('trace', 'trend', K, round(j4@teststat, digits = 3), j4@cval)
out5 <- cbind('trace', 'const', K, round(j5@teststat, digits = 3), j5@cval)
out6 <- cbind('trace', 'none',  K, round(j6@teststat, digits = 3), j6@cval)
jjci <- rbind(out1, out2, out3, out4, out5, out6)
colnames(jjci) <- c('test 1', 'test 2', 'lag', 'statistic', 
    'c.v 10%', 'c.v 5%', 'c.v 1%')
rownames(jjci) <- 1:nrow(jjci) 
(table.2 <- data.frame(jjci))

# EG cointegration
LR <- lm(formula = prVi ~ prCh); summary(LR)
(LR.coef <- round(summary(LR)$coefficients, digits = 3))
(ry <- ts(data = residuals(LR), start = start(prCh), end = end(prCh), 
  frequency = 12))
eg <- ur.df(y = ry, type = c('none'), lags = 1)
eg2 <- ur.df2(y = ry, type = c('none'), lags = 1)
(eg4 <- Box.test(eg@res, lag = 4, type = 'Ljung') )
(eg8 <- Box.test(eg@res, lag = 8, type = 'Ljung') )
(eg12 <- Box.test(eg@res, lag = 12, type = 'Ljung'))
EG.coef <- coefficients(eg@testreg)[1, 1]
EG.tval <- coefficients(eg@testreg)[1, 3]
(res.EG <- round(t(data.frame(EG.coef, EG.tval, eg2$aic, eg2$bic, 
  eg4$p.value, eg8$p.value, eg12$p.value)), digits = 3))

# -------------------------------------------------------------------------
# 4. Threshold cointegration
# best threshold 
test <- ciTarFit(y = prVi, x = prCh); test; names(test)
t3 <- ciTarThd(y = prVi, x = prCh, model = 'tar', lag = 0); plot(t3)
time.org <- proc.time() 
(th.tar <- t3$basic)
for (i in 1:12) {  # about 20 seconds
  t3a <- ciTarThd(y = prVi, x = prCh, model = 'tar', lag = i) 
  th.tar[i+2] <- t3a$basic[, 2]
}
th.tar 
time.org - proc.time()
 
t4 <- ciTarThd(y = prVi, x = prCh, model = 'mtar', lag = 0) 
(th.mtar <- t4$basic); plot(t4)
for (i in 1:12) {  # about 36 seconds
  t4a <- ciTarThd(y = prVi, x = prCh, model = 'mtar', lag = i) 
  th.mtar[i+2] <- t4a$basic[,2]
}
th.mtar

t.tar <- -8.041; t.mtar <- -0.451     # lag = 0 to 4; final choices
# t.tar <- -8.701 ; t.mtar <- -0.451  # lag = 5 to 12

mx <- 12  # lag selection
(g1 <-ciTarLag(y=prVi, x=prCh, model='tar',  maxlag = mx, thresh = 0))
(g2 <-ciTarLag(y=prVi, x=prCh, model='mtar', maxlag = mx, thresh = 0))
(g3 <-ciTarLag(y=prVi, x=prCh, model='tar',  maxlag = mx, thresh = t.tar))
(g4 <-ciTarLag(y=prVi, x=prCh, model='mtar', maxlag = mx, thresh = t.mtar))
plot(g1)

# Figure of threshold selection: mtar at lag = 3 (Figure 3 data)
(t5 <- ciTarThd(y=prVi, x=prCh, model = 'mtar', lag = 3, th.range = 0.15))
plot(t5) 

# Table 3 Results of EG and threshold cointegration combined
vv <- 3
(f1 <- ciTarFit(y=prVi, x=prCh, model = 'tar',  lag = vv, thresh = 0))
(f2 <- ciTarFit(y=prVi, x=prCh, model = 'tar',  lag = vv, thresh = t.tar ))
(f3 <- ciTarFit(y=prVi, x=prCh, model = 'mtar', lag = vv, thresh = 0))
(f4 <- ciTarFit(y=prVi, x=prCh, model = 'mtar', lag = vv, thresh = t.mtar))
 
r0 <- cbind(summary(f1)$dia, summary(f2)$dia,
            summary(f3)$dia, summary(f4)$dia)
diag <- r0[c(1:4, 6:7, 12:14, 8, 9, 11), c(1, 2, 4, 6, 8)]
rownames(diag) <- 1:nrow(diag); diag

e1 <- summary(f1)$out; e2 <- summary(f2)$out
e3 <- summary(f3)$out; e4 <- summary(f4)$out; rbind(e1, e2, e3, e4)
ee <- list(e1, e2, e3, e4); vect <- NULL
for (i in 1:4) {
  ef <- data.frame(ee[i])
  vect2 <- c(paste(ef[3, 'estimate'], ef[3, 'sign'], sep = ''),
    paste('(', ef[3, 't.value'], ')', sep = ''),
    paste(ef[4, 'estimate'], ef[4, 'sign'], sep = ''),
    paste('(', ef[4, 't.value'], ')', sep = ''))
  vect <- cbind(vect, vect2)
}
item <- c('pos.coeff','pos.t.value', 'neg.coeff','neg.t.value')
ve <- data.frame(cbind(item, vect)); colnames(ve) <- colnames(diag)
(res.CI <- rbind(diag, ve)[c(1:2, 13:16, 3:12), ])
rownames(res.CI) <- 1:nrow(res.CI)
res.CI$Engle <- '__'
res.CI[c(3, 4, 9:13), 'Engle'] <- res.EG[, 1]
res.CI[4, 6] <- paste('(', res.CI[4, 6], ')', sep = '')
(table.3 <- res.CI[, c(1, 6, 2:5)])

# -------------------------------------------------------------------------
# 5. Asymmstric error correction model
(sem <- ecmSymFit(y = prVi, x = prCh, lag = 4)); names(sem)
(aem <- ecmAsyFit(y = prVi, x = prCh, lag = 4, model = 'mtar', 
   split = TRUE, thresh = t.mtar))
(ccc <- summary(aem))
coe <- cbind(as.character(ccc[1:19, 2]), 
  paste(ccc[1:19, 'estimate'], ccc$signif[1:19], sep = ''), 
  ccc[1:19, 't.value'],
  paste(ccc[20:38, 'estimate'], ccc$signif[20:38],sep = ''), 
  ccc[20:38, 't.value']) 
colnames(coe) <- c('item', 'CH.est', 'CH.t', 'VI.est','VI.t')

(edia <- ecmDiag(aem, 3)); (ed <- edia[c(1, 6:9), ])
ed2 <- cbind(ed[, 1:2], '_', ed[, 3], '_'); colnames(ed2) <- colnames(coe)
(tes <- ecmAsyTest(aem)$out); (tes2 <- tes[c(2, 3, 5, 11:13, 1), -1])
tes3 <- cbind(as.character(tes2[, 1]), 
  paste(tes2[, 2], tes2[, 6], sep = ''), 
  paste('[', round(tes2[, 4], digits = 2), ']', sep = ''),
  paste(tes2[, 3], tes2[, 7], sep = ''), 
  paste('[', round(tes2[, 5], digits = 2), ']', sep = ''))
colnames(tes3) <- colnames(coe)
(table.4 <- data.frame(rbind(coe, ed2, tes3)))

# -------------------------------------------------------------------------
# 6. Output
(output <- listn(table.1, table.2, table.3, table.4))
write.list(z = output, file = 'OutBedTable.csv')