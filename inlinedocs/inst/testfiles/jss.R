fermat.test <- function
### Test an integer for primality using Fermat's Little Theorem.
(n ##<< The integer to test.
 ){
  a <- floor(runif(1,min=1,max=n))
  a^n %% n == a
### Whether the integer passes the Fermat test for a randomized
### \eqn{0<a<n}.
}
is.pseudoprime <- function # Check an integer for pseudo-primality.
### A number is pseudo-prime if it is probably prime, the basis of
### which is the probabalistic Fermat test; if it passes two such
### tests, the chances are better than 3 out of 4 that \eqn{n} is
### prime.
##references<< Abelson, Hal; Jerry Sussman, and Julie
##Sussman. Structure and Interpretation of Computer
##Programs. Cambridge: MIT Press, 1984.
(n,    ##<< Integer to test for pseudoprimality.
 times
### Number of Fermat tests to perform. More tests are more likely to
### give accurate results.
 ){
  if(times==0)TRUE
  ##seealso<< \code{\link{fermat.test}}
  else if(fermat.test(n)) is.pseudoprime(n,times-1)
  else FALSE
### logical TRUE if n is probably prime.
}
try.several.times <- structure(function
### Test an integer for primality using different numbers of tests.
(n,    ##<< integer to test for primality.
 times ##<< vector of number of tests to try.
 ){
  is.prime <- sapply(times,function(t)is.pseudoprime(n,t))
  ##value<< data.frame with columns:
  data.frame(times, ##<< number of Fermat tests.
             is.prime, ##<< TRUE if probably prime
             n) ##<< Integer tested.
  ##end<<
},ex=function(){
  try.several.times(6,1:5)
  try.several.times(5,1:5)
})
setClass("DocLink", # Link documentation among related functions
### The \code{DocLink} class provides the basis for hooking together
### documentation of related classes/functions/objects. The aim is that
### documentation sections missing from the child are inherited from
### the parent class.
    representation(name = "character", ##<< name of object
                   created = "character", ##<< how created
                   parent = "character", ##<< parent class or NA
                   code = "character", ##<< actual source lines
                   description = "character") ##<< preceding description block
         )

.result <- 
 list(fermat.test = list(definition = "fermat.test <- function\n### Test an integer for primality using Fermat's Little Theorem.\n(n ##<< The integer to test.\n ){\n  a <- floor(runif(1,min=1,max=n))\n  a^n %% n == a\n### Whether the integer passes the Fermat test for a randomized\n### \\eqn{0<a<n}.\n}",  
     description = "Test an integer for primality using Fermat's Little Theorem.",  
     value = "Whether the integer passes the Fermat test for a randomized\n\\eqn{0<a<n}.",  
     `item{n}` = "The integer to test.", format = "", title = "fermat test"),  
     is.pseudoprime = list(definition = "is.pseudoprime <- function # Check an integer for pseudo-primality.\n### A number is pseudo-prime if it is probably prime, the basis of\n### which is the probabalistic Fermat test; if it passes two such\n### tests, the chances are better than 3 out of 4 that \\eqn{n} is\n### prime.\n##references<< Abelson, Hal; Jerry Sussman, and Julie\n##Sussman. Structure and Interpretation of Computer\n##Programs. Cambridge: MIT Press, 1984.\n(n,    ##<< Integer to test for pseudoprimality.\n times\n### Number of Fermat tests to perform. More tests are more likely to\n### give accurate results.\n ){\n  if(times==0)TRUE\n  ##seealso<< \\code{\\link{fermat.test}}\n  else if(fermat.test(n)) is.pseudoprime(n,times-1)\n  else FALSE\n### logical TRUE if n is probably prime.\n}",  
         description = "A number is pseudo-prime if it is probably prime, the basis of\nwhich is the probabalistic Fermat test; if it passes two such\ntests, the chances are better than 3 out of 4 that \\eqn{n} is\nprime.",  
         `item{times}` = "Number of Fermat tests to perform. More tests are more likely to\ngive accurate results.",  
         value = "logical TRUE if n is probably prime.", references = "Abelson, Hal; Jerry Sussman, and Julie\nSussman. Structure and Interpretation of Computer\nPrograms. Cambridge: MIT Press, 1984.",  
         `item{n}` = "Integer to test for pseudoprimality.", seealso = "\\code{\\link{fermat.test}}",  
         title = "Check an integer for pseudo-primality.", format = ""),  
     try.several.times = list(definition = "try.several.times <- structure(function\n### Test an integer for primality using different numbers of tests.\n(n,    ##<< integer to test for primality.\n times ##<< vector of number of tests to try.\n ){\n  is.prime <- sapply(times,function(t)is.pseudoprime(n,t))\n  ##value<< data.frame with columns:\n  data.frame(times, ##<< number of Fermat tests.\n             is.prime, ##<< TRUE if probably prime\n             n) ##<< Integer tested.\n  ##end<<\n},ex=function(){\n  try.several.times(6,1:5)\n  try.several.times(5,1:5)\n})",  
         description = "Test an integer for primality using different numbers of tests.",  
         `item{n}` = "integer to test for primality.", `item{times}` = "vector of number of tests to try.",  
         value = "data.frame with columns:\n\\item{times}{number of Fermat tests.}\n\\item{is.prime}{TRUE if probably prime}\n\\item{n}{Integer tested.}",  
         format = "", title = "try several times", examples = "\ntry.several.times(6,1:5)\ntry.several.times(5,1:5)\n"),  
     `DocLink-class` = list(`item{name}` = "name of object", `item{created}` = "how created",  
         `item{parent}` = "parent class or NA", `item{code}` = "actual source lines",  
         `item{description}` = "preceding description block",  
         description = "The \\code{DocLink} class provides the basis for hooking together\ndocumentation of related classes/functions/objects. The aim is that\ndocumentation sections missing from the child are inherited from\nthe parent class.",  
         title = "Link documentation among related functions",  
         details = "Objects can be created by calls of the form \\code{new(DocLink ...)}",  
         `section{Objects from the Class}` = "Objects can be created by calls of the form \\code{new(DocLink ...)}",  
         seealso = "", alias = "DocLink", format = "")) 
