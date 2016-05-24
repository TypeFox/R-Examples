`extract` <- 
function(fun)
{
   function.list <- c()
   package.list <- c()
   variable.list <- findGlobals(fun, merge=FALSE)$variables

   loclib <- c(grep("require", as.vector(body(fun))), 
                grep("library", as.vector(body(fun))))
   for (i in loclib){
      package.list <- c(package.list, as.character(body(fun)[[i]][[2]]))
   }
   allfuns <- findGlobals(fun, merge=FALSE)$functions
   for (i in 1:length(allfuns)){
   name <- try(get(as.character(allfuns[i])), silent=TRUE)
     if(!class(name)=="try-error"){ 
        packagename <- environmentName(environment(name))
        if (packagename=="R_GlobalEnv"){
           function.list <- c(function.list, allfuns[i])
        } else{
           if (packagename!="")
           package.list <- c(package.list, packagename)
        }
   }
}
   extraction <- list(functions=function.list, packages=package.list, variables=variable.list)
   extraction
}


#    if (is.function(fun)){
#    body.fun <- body(fun)
#    } else {
#    body.fun <- fun
#    }
#       for (i in seq(along=body.fun)){
# #print(body.fun[[i]])
#          if (is.null(body.fun[[i]]) || is.na(as.character(body.fun[[i]]))) next
#          if (try(as.character(body.fun[[i]][[1]]), silent=TRUE)=="library" || 
# 	     try(as.character(body.fun[[i]][[1]]), silent=TRUE) =="require")
#          {
#             package.list <- c(package.list, as.character(body.fun[[i]][[2]]))
# 
#          } else {
# 
#            if (length(body.fun[[i]])>1){
#               extract.list <- extract(body.fun[[i]])
#               package.list <- c(package.list, extract.list$packages)
#               function.list <- c(function.list, extract.list$functions)
#            } else {
#                name <- try(get(as.character(body.fun[[i]])), silent=TRUE)
#                if(!class(name)=="try-error"){ 
#                   packagename <- environmentName(environment(name))
# 		  if (packagename=="R_GlobalEnv"){
# 		     function.list <- c(function.list, as.character(body.fun[[i]]))
# 		  } else{
#                      if (packagename!="")
#                      package.list <- c(package.list, packagename)
# 		     }
# 		  }
# 	      }
# 	   }
#       }
#    extraction <- list(functions=function.list, packages=package.list)
#    extraction
# }
# 
