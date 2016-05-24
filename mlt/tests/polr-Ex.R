
library("MASS")

mp <- polr(Sat ~ Infl, weights = Freq, data = housing)

library("mlt")

s <- as.basis(~ Infl, data = housing, remove_intercept = TRUE)
r <- as.basis(housing$Sat)
#r <- as.basis(~ Sat, data = housing, remove_intercept = TRUE,
#              contrasts.arg = list(Sat = function(n) 
#                  contr.treatment(n, base = 3)),
#              ui = diff(diag(2)), ci = 0)

m <- ctm(r, shift = s, todist = "Logi")

mod <- mlt(m, data = housing, weights = housing$Freq)

logLik(mp)
logLik(mod)

coef(mp)
mp$zeta
coef(mod)

sqrt(diag(vcov(mp)))
sqrt(diag(vcov(mod)))
