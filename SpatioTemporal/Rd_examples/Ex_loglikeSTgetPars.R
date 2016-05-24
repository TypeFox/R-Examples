##load the data
data(mesa.model)

##Compute dimensions for the data structure
dim <- loglikeSTdim(mesa.model)

##Let's create random parameter vectors ...
x <- runif( dim$nparam.cov )
names(x) <- loglikeSTnames(mesa.model, FALSE) 
x.all <- runif( dim$nparam )
names(x.all) <- loglikeSTnames(mesa.model, TRUE) 

##... and pick them apart
str( loglikeSTgetPars(x, mesa.model) )
str( loglikeSTgetPars(x.all, mesa.model) )

##Try a somewhat more interesting covariance structure
mesa.model.alt <- updateCovf(mesa.model,
                            cov.beta=list(covf=c("exp","exp2","matern"),
                              nugget=c(TRUE,FALSE,TRUE)),
                             cov.nu=list(covf="exp", nugget="type",
                               random.effect=TRUE))
##Compute dimensions for the data structure
dim <- loglikeSTdim(mesa.model.alt)

##Let's create random parameter vectors ...
x <- runif( dim$nparam.cov )
names(x) <- loglikeSTnames(mesa.model.alt, FALSE) 
x.all <- runif( dim$nparam )
names(x.all) <- loglikeSTnames(mesa.model.alt, TRUE) 

##... and pick them apart
str( loglikeSTgetPars(x, mesa.model.alt) )
str( loglikeSTgetPars(x.all, mesa.model.alt) )
