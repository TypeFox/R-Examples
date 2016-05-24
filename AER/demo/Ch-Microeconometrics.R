###################################################
### chunk number 1: setup
###################################################
options(prompt = "R> ", continue = "+  ", width = 64,
  digits = 4, show.signif.stars = FALSE, useFancyQuotes = FALSE)

options(SweaveHooks = list(onefig =   function() {par(mfrow = c(1,1))},
                           twofig =   function() {par(mfrow = c(1,2))},                           
                           threefig = function() {par(mfrow = c(1,3))},
                           fourfig =  function() {par(mfrow = c(2,2))},
			   sixfig =   function() {par(mfrow = c(3,2))}))

library("AER")

set.seed(1071)


###################################################
### chunk number 2: swisslabor-data
###################################################
data("SwissLabor")
swiss_probit <- glm(participation ~ . + I(age^2),
  data = SwissLabor, family = binomial(link = "probit"))
summary(swiss_probit)


###################################################
### chunk number 3: swisslabor-plot eval=FALSE
###################################################
## plot(participation ~ age, data = SwissLabor, ylevels = 2:1)


###################################################
### chunk number 4: swisslabor-plot-refined
###################################################
plot(participation ~ education, data = SwissLabor, ylevels = 2:1)
fm <- glm(participation ~ education + I(education^2), data = SwissLabor, family = binomial)
edu <- sort(unique(SwissLabor$education))
prop <- sapply(edu, function(x) mean(SwissLabor$education <= x))
lines(predict(fm, newdata = data.frame(education = edu), type = "response") ~ prop, col = 2)

plot(participation ~ age, data = SwissLabor, ylevels = 2:1)
fm <- glm(participation ~ age + I(age^2), data = SwissLabor, family = binomial)
ag <- sort(unique(SwissLabor$age))
prop <- sapply(ag, function(x) mean(SwissLabor$age <= x))
lines(predict(fm, newdata = data.frame(age = ag), type = "response") ~ prop, col = 2)


###################################################
### chunk number 5: effects1
###################################################
fav <- mean(dnorm(predict(swiss_probit, type = "link")))
fav * coef(swiss_probit)


###################################################
### chunk number 6: effects2
###################################################
av <- colMeans(SwissLabor[, -c(1, 7)])
av <- data.frame(rbind(swiss = av, foreign = av),
  foreign = factor(c("no", "yes")))
av <- predict(swiss_probit, newdata = av, type = "link")
av <- dnorm(av)
av["swiss"] * coef(swiss_probit)[-7]


###################################################
### chunk number 7: effects3
###################################################
av["foreign"] * coef(swiss_probit)[-7]


###################################################
### chunk number 8: mcfadden
###################################################
swiss_probit0 <- update(swiss_probit, formula = . ~ 1)
1 - as.vector(logLik(swiss_probit)/logLik(swiss_probit0))


###################################################
### chunk number 9: confusion-matrix
###################################################
table(true = SwissLabor$participation,
  pred = round(fitted(swiss_probit)))


###################################################
### chunk number 10: confusion-matrix1
###################################################
tab <- table(true = SwissLabor$participation,
  pred = round(fitted(swiss_probit)))
tabp <- round(100 * c(tab[1,1] + tab[2,2], tab[2,1] + tab[1,2])/sum(tab), digits = 2)


###################################################
### chunk number 11: roc-plot eval=FALSE
###################################################
## library("ROCR")
## pred <- prediction(fitted(swiss_probit),
##   SwissLabor$participation)
## plot(performance(pred, "acc"))
## plot(performance(pred, "tpr", "fpr"))
## abline(0, 1, lty = 2)


###################################################
### chunk number 12: roc-plot1
###################################################
library("ROCR")
pred <- prediction(fitted(swiss_probit),
  SwissLabor$participation)
plot(performance(pred, "acc"))
plot(performance(pred, "tpr", "fpr"))
abline(0, 1, lty = 2)


###################################################
### chunk number 13: rss
###################################################
deviance(swiss_probit)
sum(residuals(swiss_probit, type = "deviance")^2)
sum(residuals(swiss_probit, type = "pearson")^2)


###################################################
### chunk number 14: coeftest eval=FALSE
###################################################
## coeftest(swiss_probit, vcov = sandwich)


###################################################
### chunk number 15: murder
###################################################
data("MurderRates")
murder_logit <- glm(I(executions > 0) ~ time + income +
  noncauc + lfp + southern, data = MurderRates,
  family = binomial)


###################################################
### chunk number 16: murder-coeftest
###################################################
coeftest(murder_logit)


###################################################
### chunk number 17: murder2
###################################################
murder_logit2 <- glm(I(executions > 0) ~ time + income +
  noncauc + lfp + southern, data = MurderRates,
  family = binomial, control = list(epsilon = 1e-15,
  maxit = 50, trace = FALSE))


###################################################
### chunk number 18: murder2-coeftest
###################################################
coeftest(murder_logit2)


###################################################
### chunk number 19: separation
###################################################
table(I(MurderRates$executions > 0), MurderRates$southern)


###################################################
### chunk number 20: countreg-pois
###################################################
data("RecreationDemand")
rd_pois <- glm(trips ~ ., data = RecreationDemand,
  family = poisson)


###################################################
### chunk number 21: countreg-pois-coeftest
###################################################
coeftest(rd_pois)


###################################################
### chunk number 22: countreg-pois-logLik
###################################################
logLik(rd_pois)


###################################################
### chunk number 23: countreg-odtest1
###################################################
dispersiontest(rd_pois)


###################################################
### chunk number 24: countreg-odtest2
###################################################
dispersiontest(rd_pois, trafo = 2)


###################################################
### chunk number 25: countreg-nbin
###################################################
library("MASS")
rd_nb <- glm.nb(trips ~ ., data = RecreationDemand)
coeftest(rd_nb)
logLik(rd_nb)


###################################################
### chunk number 26: countreg-se
###################################################
round(sqrt(rbind(diag(vcov(rd_pois)),
  diag(sandwich(rd_pois)))), digits = 3)


###################################################
### chunk number 27: countreg-sandwich
###################################################
coeftest(rd_pois, vcov = sandwich)


###################################################
### chunk number 28: countreg-OPG
###################################################
round(sqrt(diag(vcovOPG(rd_pois))), 3)


###################################################
### chunk number 29: countreg-plot
###################################################
plot(table(RecreationDemand$trips), ylab = "")


###################################################
### chunk number 30: countreg-zeros
###################################################
rbind(obs = table(RecreationDemand$trips)[1:10], exp = round(
  sapply(0:9, function(x) sum(dpois(x, fitted(rd_pois))))))


###################################################
### chunk number 31: countreg-pscl
###################################################
library("pscl")


###################################################
### chunk number 32: countreg-zinb
###################################################
rd_zinb <- zeroinfl(trips ~ . | quality + income,
  data = RecreationDemand, dist = "negbin")


###################################################
### chunk number 33: countreg-zinb-summary
###################################################
summary(rd_zinb)


###################################################
### chunk number 34: countreg-zinb-expected
###################################################
round(colSums(predict(rd_zinb, type = "prob")[,1:10]))


###################################################
### chunk number 35: countreg-hurdle
###################################################
rd_hurdle <- hurdle(trips ~ . | quality + income,
  data = RecreationDemand, dist = "negbin")
summary(rd_hurdle)


###################################################
### chunk number 36: countreg-hurdle-expected
###################################################
round(colSums(predict(rd_hurdle, type = "prob")[,1:10]))


###################################################
### chunk number 37: tobit1
###################################################
data("Affairs")
aff_tob <- tobit(affairs ~ age + yearsmarried +
  religiousness + occupation + rating, data = Affairs)
summary(aff_tob)


###################################################
### chunk number 38: tobit2
###################################################
aff_tob2 <- update(aff_tob, right = 4)
summary(aff_tob2)


###################################################
### chunk number 39: tobit3
###################################################
linearHypothesis(aff_tob, c("age = 0", "occupation = 0"),
  vcov = sandwich)


###################################################
### chunk number 40: numeric-response
###################################################
SwissLabor$partnum <- as.numeric(SwissLabor$participation) - 1


###################################################
### chunk number 41: kleinspady eval=FALSE
###################################################
## library("np")
## swiss_bw <- npindexbw(partnum ~ income + age + education +
##   youngkids + oldkids + foreign + I(age^2), data = SwissLabor,
##   method = "kleinspady", nmulti = 5)


###################################################
### chunk number 42: kleinspady-bw eval=FALSE
###################################################
## summary(swiss_bw)


###################################################
### chunk number 43: kleinspady-summary eval=FALSE
###################################################
## swiss_ks <- npindex(bws = swiss_bw, gradients = TRUE)
## summary(swiss_ks)


###################################################
### chunk number 44: probit-confusion
###################################################
table(Actual = SwissLabor$participation, Predicted = 
  round(predict(swiss_probit, type = "response")))


###################################################
### chunk number 45: bw-tab
###################################################
data("BankWages")
edcat <- factor(BankWages$education)
levels(edcat)[3:10] <- rep(c("14-15", "16-18", "19-21"),
  c(2, 3, 3))
tab <- xtabs(~ edcat + job, data = BankWages)
prop.table(tab, 1)


###################################################
### chunk number 46: bw-plot eval=FALSE
###################################################
## plot(job ~ edcat, data = BankWages, off = 0)


###################################################
### chunk number 47: bw-plot1
###################################################
plot(job ~ edcat, data = BankWages, off = 0)
box()


###################################################
### chunk number 48: bw-multinom
###################################################
library("nnet")
bank_mnl <- multinom(job ~ education + minority,
  data = BankWages, subset = gender == "male", trace = FALSE)


###################################################
### chunk number 49: bw-multinom-coeftest
###################################################
coeftest(bank_mnl)


###################################################
### chunk number 50: bw-polr
###################################################
library("MASS")
bank_polr <- polr(job ~ education + minority, 
  data = BankWages, subset = gender == "male", Hess = TRUE)
coeftest(bank_polr)


###################################################
### chunk number 51: bw-AIC
###################################################
AIC(bank_mnl)
AIC(bank_polr) 


