
# Forces evaluation of objects with names in "toForce" from "envir" enviroment to "cluster".
# If ALL == TRUE it evaluates all the objects in "envir"

.forceEval <- function(toForce = c(), envir = parent.frame(), ALL = FALSE)
{
  allNames <- toForce
  
  if(ALL) allNames <- c(allNames, ls(envir = envir))
  
  lapply(allNames, function(input) force(get(input, envir = envir)) )
  
  return( NULL )
}


# Test

# f <- function(x, y){ function() c(x, y) }
# lf <- vector("list", 5)
# for (i in seq_along(lf)) lf[[i]] <- f(i, i)
# lf[[1]]()  # returns 5, 5 which is wrong
# 
# f <- function(x, y){ .forceEval(c("x", "y")); function() c(x, y) }
# lf <- vector("list", 5)
# for (i in seq_along(lf)) lf[[i]] <- f(i, i)
# lf[[1]]()  # returns 1, 1 which right