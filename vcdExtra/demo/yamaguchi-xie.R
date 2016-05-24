## Models for Yamaguchi (1987) data on social mobility in US, UK and Japan, following Xie (1992)
## These models reproduce the results in Table 1, appplied to the off-diagonal cells

library(gnm)
library(vcdExtra)

data(Yamaguchi87)
# create table form
Yama.tab <- xtabs(Freq ~ Father + Son + Country, data=Yamaguchi87)

# define labeling_args for convenient reuse in 3-way displays
largs <- list(rot_labels=c(right=0), offset_varnames = c(right = 0.6), offset_labels = c(right = 0.2),
		set_varnames = c(Father="Father's status", Son="Son's status") 
)

# no association between F and S given country ('perfect mobility')
# asserts same associations for all countries

yamaNull <- gnm(Freq ~ (Father + Son) * Country, data=Yamaguchi87, family=poisson)
LRstats(yamaNull)
mosaic(yamaNull, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs,
		main="[FC][SC] Null [FS] association (perfect mobility)")

## same, with data in xtabs form
#yamaNull <- gnm(Freq ~ (Father + Son) * Country, data=Yama.tab, family=poisson)
#LRstats(yamaNull)
#mosaic(yamaNull, ~Country + Son + Father, condvars="Country", 
#	labeling_args=largs,
#	main="[FC][SC] Null [FS] association (perfect mobility)")

# ignore diagonal cells, overall
#yamaDiag0 <- gnm(Freq ~ (Father + Son) * Country + Diag(Father, Son), data=Yama.tab, family=poisson)
#LRstats(yamaDiag0)
# same, using update()
yamaDiag0 <- update(yamaNull, ~ . + Diag(Father, Son))
LRstats(yamaDiag0)

# ignore diagonal cells in each Country [Model NA in Xie(1992), Table 1]
yamaDiag <- update(yamaNull, ~ . + Diag(Father, Son):Country)
LRstats(yamaDiag)
mosaic(yamaDiag, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="[FC][SC] Quasi perfect mobility, +Diag(F,S)")


# fit models using integer scores for rows/cols 
Rscore <- as.numeric(Yamaguchi87$Father)
Cscore <- as.numeric(Yamaguchi87$Son)

# cross-nationally homogeneous row effect associations (Xie, model R_o)
yamaRo <- update(yamaDiag, ~ . + Father:Cscore)
LRstats(yamaRo)
mosaic(yamaRo, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model Ro: homogeneous row effects, +Father:j ")

# cross-nationally log multiplicative row effect associations (Xie, model R_x)
yamaRx <- update(yamaDiag, ~ . + Mult(Father:Cscore, Exp(Country)))
LRstats(yamaRx)

# cross-nationally homogeneous col effect associations (Xie, model C_o)
yamaCo <- update(yamaDiag, ~ . + Rscore:Son)
LRstats(yamaCo)
mosaic(yamaCo, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model Co: homogeneous col effects, +i:Son")


# cross-nationally log multiplicative col effect associations (Xie, model C_x)
yamaCx <- update(yamaDiag, ~ . + Mult(Rscore:Son, Exp(Country)))
LRstats(yamaCx)

# cross-nationally homogeneous row + col effect associations I (Xie, model (R+C)_o)
yamaRpCo <- update(yamaDiag, ~ . + Father:Cscore + Rscore:Son)
LRstats(yamaRpCo)
mosaic(yamaRpCo, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model (R+C)o: homogeneous, F:j + i:S")

# cross-nationally log multiplicative row + col effect associations I (Xie, model (R+C)_x)
yamaRpCx <- update(yamaDiag, ~ . + Mult(Father:Cscore + Rscore:Son, Exp(Country)))
LRstats(yamaRpCx)
mosaic(yamaRpCx, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model (R+C)x: log multiplicative (Fj + iS) : Country")

# cross-nationally homogeneous row and col effect associations II (Xie, model RC_o)
yamaRCo <- update(yamaDiag, ~ . + Mult(Father,Son))
LRstats(yamaRCo)
mosaic(yamaRCo, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model RCo: homogeneous RC(1)")

# cross-nationally log multiplicative row and col effect associations II (Xie, model RC_x)
yamaRCx <- update(yamaDiag, ~ . + Mult(Father,Son, Exp(Country)))
LRstats(yamaRCx)
mosaic(yamaRCx, ~Country + Son + Father, condvars="Country", 
		labeling_args=largs, gp=shading_Friendly,
		main="Model RCx: log multiplicative RC(1) : Country")

# cross-nationally homogeneous full two-way RxC association (Xie, model FI_o)
yamaFIo <- update(yamaDiag, ~ . + Father:Son)
LRstats(yamaFIo)

# cross-nationally log multiplicative full two-way RxC association (Xie, model FI_x)
yamaFIx <- update(yamaDiag, ~ . + Mult(Father:Son, Exp(Country)))
LRstats(yamaFIx)

# compare models
models <- glmlist(yamaNull, yamaDiag, yamaRo, yamaRx, yamaCo, yamaCx, yamaRpCo, yamaRpCx, yamaRCo, yamaRCx, yamaFIo, yamaFIx)
LRstats(models)


# extract models sumaries, consider as factorial of RC model by layer model

BIC <- matrix(LRstats(models)$BIC[-(1:2)], 5, 2, byrow=TRUE)
dimnames(BIC) <- list("Father-Son model" = c("row eff", "col eff", "row+col", "RC(1)", "R:C"),
		"Country model" = c("homogeneous", "log multiplicative"))
BIC

matplot(BIC, type='b', xlab="Father-Son model", xaxt='n', pch=15:16, cex=1.5, cex.lab=1.5,
		main="Yamaguchi-Xie models: R:C model by Layer model Summary")
axis(side=1, at=1:nrow(BIC), labels=rownames(BIC), cex.axis=1.2)
text(5, BIC[5,], colnames(BIC), pos=2, col=1:2, cex=1.2)
text(5, max(BIC[5,])+10, "Country model", pos=2, cex=1.3)

AIC <- matrix(LRstats(models)$AIC[-(1:2)], 5, 2, byrow=TRUE)
dimnames(AIC) <- list("Father-Son model" = c("row eff", "col eff", "row+col", "RC(1)", "R:C"),
		"Country model" = c("homogeneous", "log multiplicative"))
AIC

matplot(AIC, type='b', xlab="Father-Son model", xaxt='n', pch=15:16, cex=1.5, cex.lab=1.5,
		main="Yamaguchi-Xie models: R:C model by Layer model Summary")
axis(side=1, at=1:nrow(AIC), labels=rownames(AIC), cex.axis=1.2)
text(5, AIC[5,], colnames(AIC), pos=2, col=1:2, cex=1.2)
text(5, max(AIC[5,])+10, "Country model", pos=2, cex=1.3)

