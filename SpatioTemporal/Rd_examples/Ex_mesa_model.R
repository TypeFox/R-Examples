##load the data
data(mesa.model)

##examine components
names(mesa.model)
print(mesa.model)
summary(mesa.model)

##requested geographic and spatio-temporal covariates
mesa.model$LUR.list
mesa.model$ST.list

##covariates for the temporal intercept
head(mesa.model$LUR$const)
##...and the two smooth temporal trends
head(mesa.model$LUR$V1)
head(mesa.model$LUR$V2)

##Some important dimensions of the model
loglikeSTdim(mesa.model)
