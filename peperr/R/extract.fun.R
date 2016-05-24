`extract.fun` <-
function(funs=NULL) 
{
   require(codetools)
   if (is.null(funs)){
      extract.packages <- (.packages())
      variables <- c()
      for (i in seq(along=ls(name=".GlobalEnv"))){
          if (!is.function(get(ls(name=".GlobalEnv")[i]))){
             variables <- c(variables, ls(name=".GlobalEnv")[i])
          }
      }
      functions <- c()
      for (i in seq(along=ls(name=".GlobalEnv"))){
         name <- try(get(as.character(ls(name=".GlobalEnv")[i])), silent=TRUE)
         packagename <- environmentName(environment(name))
         if (packagename=="R_GlobalEnv"){
            functions <- c(functions, ls(name=".GlobalEnv")[i])
         }
      }
   } else {
      packages <- c()
      functions <- c()
      variables <- c()
      for (j in seq(along=funs)){
         if (is.function(funs[[j]])){
            packages <- c(packages, extract(funs[[j]])$packages)
            functions <- c(functions, extract(funs[[j]])$functions) 
            variables <- c(variables, extract(funs[[j]])$variables)
         }
      }
   extract.packages <- unique(packages)
   }

   if (length(extract.packages)==0) extract.packages <- NULL

   extract.functions <- unique(functions)
   extract.variables <- unique(variables)
    
   extract.fun.list <- list(packages=extract.packages, 
        functions=extract.functions, variables=extract.variables)
   extract.fun.list
}

