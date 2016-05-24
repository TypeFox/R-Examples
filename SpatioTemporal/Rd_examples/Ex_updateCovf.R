##load the data
data(mesa.model)

##covariance specification:
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=TRUE, random.effect=FALSE)

##Simple covariance structure
updateCovf(mesa.model, cov.beta, cov.nu)

##different behaviour for different beta:s
cov.beta <- list(covf=c("exp","exp2","matern"), nugget=c(TRUE,FALSE,FALSE))
updateCovf(mesa.model, cov.beta, cov.nu)

##Spatially varying nugget
cov.nu <- list(covf="exp", nugget="type", random.effect=FALSE)
print(tmp <- updateCovf(mesa.model, cov.beta, cov.nu))
##lets study the regression matrix for the nugget
str(tmp$cov.nu$nugget.matrix)
head(tmp$cov.nu$nugget.matrix)
