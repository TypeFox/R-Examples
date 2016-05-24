"lm.beta" <-
function (MOD) 
{
b <- summary(MOD)$coef[-1,1]
sx <- sapply(MOD$model[-1], sd)
sy <- sapply(MOD$model[1], sd)
beta <- b * sx /  sy
return(beta)
}

