# A. Without S3 or S4
dog.a <- function(w) {
  if (!is.numeric(w)) stop("w has to be numeric.\n")
  z <- w / 100
  m <- ifelse(test = w[1] < 50, yes = "small", no = "big")
  result <- list(w = w, z = z, m = m)
  return(result)
}
ra <- dog.a(w = 1:8); str(ra)
print(ra)  # print.default invoked; display everything
plot(ra)   # plot.default does not work for this example

# Manually without S3/S4: selected printing and plotting
ra$z
plot(x = ra$w, y = ra$z, col = "purple", main = "Drawing a graph manually")

# B. Object-oriented programming with S3 ----------------------------------
dog.b <- function(w) {
  if (!is.numeric(w)) stop("w has to be numeric.\n")
  z <- w / 100
  m <- ifelse(test = w[1] < 50, yes = "small", no = "big")
  result <- list(w = w, z = z, m = m)
  class(result) <- "ego"  # S3 class assigned
  return(result)
}

print.ego <- function(x, ...) {  # define an S3 print method
  cat("\n=== My personalized output by S3 =========\n")
  print(x$z)
  invisible(x)
}

plot.ego <- function(x, ...) {  # define an S3 plot method
  plot(x = x$w, y = x$z, col ="blue", main = "My first graph by S3!")
  invisible(NULL)  
}

rb <- dog.b(w = 1:8); str(rb)
print(rb)  # S3 print.ego invoked; selected printing and formatting
plot(rb)   # S3 plot.ego invoked  

# C. Object-oriented programming with S4 ----------------------------------
setClass(Class = 'earth',  # define a new S4 class 
  slots = c(w = 'numeric', z = 'numeric', m = 'character')
)

dog.c <- function(w) {
  if (!is.numeric(w)) stop("w has to be numeric.\n")
  z <- w / 100
  m <- ifelse(test = w[1] < 50, yes = "small", no = "big")
  result <- new("earth", w = w, z = z, m = m)  # create an S4 object
  return(result)
}

setMethod(f = "show",  # new S4 print method  
  signature = signature(object = "earth"),  
  definition = function(object) {
    cat('\n=== My personalized output by S4 =========\n')
    print(object@z)
    invisible(object)
  }
)

setMethod(f = "plot",  # new S4 plot method 
  signature = signature(x = "earth"),  
  definition = function(x) {
    plot(x = x@w, y = x@z, col = "red", main = "My first graph by S4!")
    invisible(NULL)
  }  
)

rc <- dog.c(w = 1:8); str(rc)
show(rc)  # S4 show method for 'earth' invoked; selected print + format
plot(rc)  # S4 plot method for 'earth' invoked