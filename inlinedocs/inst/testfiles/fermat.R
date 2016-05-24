fermat.test <- function#Test an integer for primality with Fermat's little theorem.
### Fermat's little theorem states that if \eqn{n} is a prime number
### and \eqn{a} is any positive integer less than \eqn{n}, then
### \eqn{a} raised to the \eqn{n}th power is congruent to \eqn{a\
### modulo\ n}{a modulo n}.
##references<< \url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}
(n ##<< the integer to test for primality.
 ){
  a <- floor(runif(1,min=1,max=n))
  ##note<< \code{fermat.test} doesn't work for integers above
  ##approximately 15 because modulus loses precision.
  a^n %% n == a
### Whether the integer passes the Fermat test for a randomized
### \eqn{0<a<n}
}

is.pseudoprime <- function#Check an integer for pseudo-primality to an arbitrary precision.
### A number is pseudo-prime if it is probably prime, the basis of
### which is the probabalistic Fermat test; if it passes two such
### tests, the chances are better than 3 out of 4 that \eqn{n} is
### prime.
##references<< Abelson, Hal; Jerry Sussman, and Julie
##Sussman. Structure and Interpretation of Computer
##Programs. Cambridge: MIT Press, 1984.
(n, ##<< the integer to test for pseudoprimality.
 times ##<< the number of Fermat tests to perform
){
  if(times==0)TRUE
  ##seealso<< \code{\link{fermat.test}}
  else if(fermat.test(n)) is.pseudoprime(n,times-1)
  else FALSE
### Whether the number is pseudoprime
}

.result <- list(fermat.test=list(format="",
  definition="fermat.test <- function#Test an integer for primality with Fermat's little theorem.\n### Fermat's little theorem states that if \\eqn{n} is a prime number\n### and \\eqn{a} is any positive integer less than \\eqn{n}, then\n### \\eqn{a} raised to the \\eqn{n}th power is congruent to \\eqn{a\\\n### modulo\\ n}{a modulo n}.\n##references<< \\url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}\n(n ##<< the integer to test for primality.\n ){\n  a <- floor(runif(1,min=1,max=n))\n  ##note<< \\code{fermat.test} doesn't work for integers above\n  ##approximately 15 because modulus loses precision.\n  a^n %% n == a\n### Whether the integer passes the Fermat test for a randomized\n### \\eqn{0<a<n}\n}",
  title="Test an integer for primality with Fermat's little theorem.",
  description="Fermat's little theorem states that if \\eqn{n} is a prime number\nand \\eqn{a} is any positive integer less than \\eqn{n}, then\n\\eqn{a} raised to the \\eqn{n}th power is congruent to \\eqn{a\\\nmodulo\\ n}{a modulo n}.",
  value="Whether the integer passes the Fermat test for a randomized\n\\eqn{0<a<n}",
  references="\\url{http://en.wikipedia.org/wiki/Fermat's_little_theorem}",
  `item{n}`="the integer to test for primality.",
  note="\\code{fermat.test} doesn't work for integers above\napproximately 15 because modulus loses precision."),
                is.pseudoprime=list(format="",
  definition="is.pseudoprime <- function#Check an integer for pseudo-primality to an arbitrary precision.\n### A number is pseudo-prime if it is probably prime, the basis of\n### which is the probabalistic Fermat test; if it passes two such\n### tests, the chances are better than 3 out of 4 that \\eqn{n} is\n### prime.\n##references<< Abelson, Hal; Jerry Sussman, and Julie\n##Sussman. Structure and Interpretation of Computer\n##Programs. Cambridge: MIT Press, 1984.\n(n, ##<< the integer to test for pseudoprimality.\n times ##<< the number of Fermat tests to perform\n){\n  if(times==0)TRUE\n  ##seealso<< \\code{\\link{fermat.test}}\n  else if(fermat.test(n)) is.pseudoprime(n,times-1)\n  else FALSE\n### Whether the number is pseudoprime\n}",
  title="Check an integer for pseudo-primality to an arbitrary precision.",
  description="A number is pseudo-prime if it is probably prime, the basis of\nwhich is the probabalistic Fermat test; if it passes two such\ntests, the chances are better than 3 out of 4 that \\eqn{n} is\nprime.",
  `item{times}`= "the number of Fermat tests to perform",
  value="Whether the number is pseudoprime",
  references="Abelson, Hal; Jerry Sussman, and Julie\nSussman. Structure and Interpretation of Computer\nPrograms. Cambridge: MIT Press, 1984.",
  `item{n}`="the integer to test for pseudoprimality.",
  seealso="\\code{\\link{fermat.test}}"))
