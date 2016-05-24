if(length(find.package(c("survival", "coxme"), quiet = TRUE)) == 2) {

library(coxme)
library(MuMIn)
options(na.action = "na.fail")


lung$temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))
lung1 <- na.omit(lung)

rfit0 <- coxme(Surv(time, status) ~ ph.ecog * ph.karno + (age | 1) + (wt.loss | 1),
	data = lung1)

getAllTerms(rfit0)
	
stopifnot(is.numeric(coeffs(rfit0)))

dd <- dredge(rfit0, eval = TRUE, trace = TRUE)
coeffs(dd)
model.sel(dd, rank = AIC)
summary(ma <- model.avg(dd))

library(survival)

fm0 <- coxph(formula = Surv(time, status) ~ 1, data = lung1)
fm1 <- coxph(formula = Surv(time, status) ~ temp, data = lung1)
fme <- coxme(formula = Surv(time, status) ~ ph.ecog + (temp | 1), data = lung1)
fme1 <- coxme(formula = Surv(time, status) ~ ph.ecog + (age | 1) + (wt.loss | 1), data = lung1)
fme2 <- coxme(formula = Surv(time, status) ~ ph.ecog + (wt.loss | 1), data = lung1)
fm2 <- coxph(formula = Surv(time, status) ~ temp *  ph.ecog, data = lung1)

avg <- model.avg(dredge(fm2), fit = T)

#avgpred(avg, type = "response")


(ms <- model.sel(fm0, fm1, fme, fme1, fme2))

summary(avg <- model.avg(ms))

predict(fm1, se.fit = T)


if(length(find.package("nlme", quiet = TRUE)) == 1) {
	#fit1 <- lme(effort ~ Type, data=ergoStool, random= ~1|Subject/ran1, method="ML")

	if(exists("lmekin", mode = "function", envir = asNamespace("coxme"))) {
		data("ergoStool", package = "nlme")
		fit2 <- lmekin(effort ~ Type + (1|Subject), data = ergoStool)
		dd <- dredge(fit2, trace = TRUE)
		summary(ma <- model.avg(dd))
	}
}

}