setOldClass("family")
setOldClass("mcmc.list")
setOldClass("polr")
setOldClass("bugs")
setOldClass("svyglm")

setClass("balance",
     representation(
            rawdata = "data.frame",
            matched = "data.frame",
            factor = "logical")
)




setClass("bayesglm",
     representation(
            formula = "formula",
            family = "family",
            prior.mean = "numeric", 
            prior.scale = "numeric", 
            prior.df = "numeric"),
    contains = "glm"
)


#setClass("bayesglm.h",
#     representation(
#            formula = "formula",
#            family = "family",
#            prior.mean = "numeric", 
#            prior.scale = "numeric", 
#            prior.df = "numeric",
#            batch = "numeric"),
#    contains = "bayesglm"
#)

#setClass("polr",
#     representation(
#            formula = "formula",
#            Hess = "logical",
#            method = "character"
##            prior.mean = "numeric", 
##            prior.scale = "numeric", 
##            prior.df = "numeric",
##            prior.mean.for.cutpoints = "numeric", 
##            prior.scale.for.cutpoints = "numeric",  
##            prior.df.for.cutpoints = "numeric"
#            ),
#    contains="oldClass"
#)


setClass("bayespolr",
     representation(
            formula = "formula",
            Hess = "logical",
            method = "character",
            prior.mean = "numeric", 
            prior.scale = "numeric", 
            prior.df = "numeric",
            prior.mean.for.cutpoints = "numeric", 
            prior.scale.for.cutpoints = "numeric",  
            prior.df.for.cutpoints = "numeric"),
     contains = "polr"
)


setClass("sim",
     representation(
            coef = "matrix",
            sigma = "numeric")
)

setClass("sim.polr",
     representation(
            coef = "matrix",
            zeta = "matrix")
)



setClass("sim.merMod",
     representation(
            fixef = "matrix",
            ranef = "list",
            sigma = "ANY")
)


setClass("GO")
