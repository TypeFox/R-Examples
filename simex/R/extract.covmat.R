extract.covmat <-
function (model) {
type <- class(model)[1]
sum.model <- summary(model)
switch(type,
"glm" = covmat <- sum.model$cov.scaled,
"lm"  = covmat <- sum.model$cov.unscaled * sum.model$sigma^2,
"gam" = covmat <- model$Vp,
"nls" = covmat <- sum.model$cov.unscaled * sum.model$sigma^2,
"lme" = covmat <- model$apVar,
"nlme"= covmat <- model$apVar
)
return(covmat)
}

