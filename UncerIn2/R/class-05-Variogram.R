## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian

setClass(Class = "Vario", 
         representation=representation(
           model = "character",
           psill = "numeric",
           range = "numeric",
           kappa = "numeric",
           ang1 = "numeric",
           ang2 = "numeric",
           ang3 = "numeric",
           anis1 = "numeric",
           anis2 = "numeric"
         ),
         
         prototype=prototype(
           model = NA_character_,
           psill = c(),
           range = c(),
           kappa = c(),
           ang1 = c(),
           ang2 = c(),
           ang3 = c(),
           anis1 = c(),
           anis2 = c()
         )
)


setClass(Class = "Variogram", 
         representation=representation(
           minimal = "Vario",
           modalValue = "Vario",
           maximal = "Vario"
         ),
         
)


setMethod("show",
          signature(object = "Variogram"),
          definition = function(object)           
          {
            
            cat("Miminimal Variogram \n")
            cat("model \t psill \t range \t kappa \n")
            cat(paste(object@minimal@model[1],object@minimal@psill[1],object@minimal@range[1],object@minimal@kappa[1],"\n", sep="\t"))
            cat(paste(object@minimal@model[2],object@minimal@psill[2],object@minimal@range[2],object@minimal@kappa[2],"\n", sep="\t"))
            cat("-----------------------------------------\n")
            
            cat("Modal Variogram \n")
            cat("model \t psill \t range \t kappa \n")
            cat(paste(object@modalValue@model[1],object@modalValue@psill[1],object@modalValue@range[1],object@modalValue@kappa[1],"\n", sep="\t"))
            cat(paste(object@modalValue@model[2],object@modalValue@psill[2],object@modalValue@range[2],object@modalValue@kappa[2],"\n", sep="\t"))
            cat("-----------------------------------------\n")
            
            cat("Maximal Variogram  \n")
            cat("model \t psill \t range \t kappa \n")
            cat(paste(object@maximal@model[1],object@maximal@psill[1],object@maximal@range[1],object@maximal@kappa[1],"\n", sep="\t"))
            cat(paste(object@maximal@model[2],object@maximal@psill[2],object@maximal@range[2],object@maximal@kappa[2],"\n", sep="\t"))
            cat("-----------------------------------------\n")
          }
)