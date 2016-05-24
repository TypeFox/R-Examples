## ---- results='asis', echo=FALSE-----------------------------------------
cat(gsub("\\n   ", "", packageDescription("modules", fields = "Description")))

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("modules")

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("wahani/modules")

## ------------------------------------------------------------------------
library(modules)
m <- module({
  boringFunction <- function() cat("boring output")
})

m$boringFunction()

## ----error=TRUE----------------------------------------------------------
hey <- "hey"
m <- module({
  isolatedFunction <- function() hey
})
m$isolatedFunction()

## ------------------------------------------------------------------------
m <- module({
  functionWithDep <- function(x) stats::median(x)
})
m$functionWithDep(1:10)

## ------------------------------------------------------------------------
m <- module({
 
  import(stats, median) # make median from package stats available
  
  functionWithDep <- function(x) median(x)

})
m$functionWithDep(1:10)

## ------------------------------------------------------------------------
m <- module({
  
  import(stats)
  
  functionWithDep <- function(x) median(x)

})
m$functionWithDep(1:10)

## ------------------------------------------------------------------------
m <- module({
  
  export("fun")
  
  .privateFunction <- identity
  privateFunction <- identity
  fun <- identity
  
})

names(m)

## ----error=TRUE----------------------------------------------------------
library(parallel)
dependency <- identity
fun <- function(x) dependency(x) 

cl <- makeCluster(2)
clusterMap(cl, fun, 1:2)
stopCluster(cl)

## ------------------------------------------------------------------------
m <- module({
  dependency <- identity
  fun <- function(x) dependency(x) 
})

cl <- makeCluster(2)
clusterMap(cl, m$fun, 1:2)
stopCluster(cl)

## ------------------------------------------------------------------------
code <- "
import(methods)
import(aoos)
# This is an example where we rely on functions in 'aoos':
list : generic(x) %g% standardGeneric('generic')
generic(x ~ ANY) %m% as.list(x)
"

fileName <- tempfile(fileext = ".R")
writeLines(code, fileName)

## ------------------------------------------------------------------------
someModule <- use(fileName)
someModule$generic(1:2)

## ------------------------------------------------------------------------
m <- module({
  fun <- function(x) {
    ## A function for illustrating documentation
    ## x (numeric)
    x
  }
})

## ------------------------------------------------------------------------
m
m$fun

## ------------------------------------------------------------------------
mutableModule <- module({
  .num <- NULL
  get <- function() .num
  set <- function(val) .num <<- val
})
mutableModule$get()
mutableModule$set(2)

## ------------------------------------------------------------------------
complectModule <- module({
  use(.GlobalEnv$mutableModule, attach = TRUE)
  getNum <- function() get()
  set(3)
})
mutableModule$get()
complectModule$getNum()

## ------------------------------------------------------------------------
complectModule <- module({
  use(.GlobalEnv$mutableModule, attach = TRUE, reInit = FALSE)
  getNum <- function() get()
  set(3)
})
mutableModule$get()
complectModule$getNum()

## ------------------------------------------------------------------------
complectModule <- module({
  expose(.GlobalEnv$mutableModule, reInit = TRUE)
  set(4)
})
mutableModule$get()
complectModule$get()

## ------------------------------------------------------------------------
complectModule <- module({
  expose(.GlobalEnv$mutableModule, reInit = FALSE)
  set(1)
})
mutableModule$get()
complectModule$get()

## ------------------------------------------------------------------------
m <- module({
  
  import("stats", "median")
  import("modules", "module")
  
  anotherModule <- module({
    fun <- function(x) median(x)
  })
  
})

m$anotherModule$fun(1:2)

## ------------------------------------------------------------------------
m <- module({
  generic <- function(x) UseMethod("generic")
  generic.numeric <- function(x) cat("method for x ~ numeric")
})
# m$generic(1) # this won't work
use(m, attach = TRUE)
m$generic(1)

## ------------------------------------------------------------------------
m <- module({
  import("methods")
  import("aoos")
  gen(x) %g% cat("default method")
  gen(x ~ numeric) %m% cat("method for x ~ numeric")
})
m$gen("Hej")
m$gen(1)

## ------------------------------------------------------------------------
m <- module({
  import("methods")
  import("aoos")
  numeric : NewType() %type% .Object
})
m$NewType(1)

