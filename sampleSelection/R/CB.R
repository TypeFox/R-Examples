## C(B) = -(pnorm(B)*dnorm(B)*B + dnorm(B)^2)/pnorm(B)
##
## Numerically robust method recommended by Dimitrios Rizopoulos,
## KULeuven.  How to do it better?  How to prove the limit value?
CB <- function(x) {
   ifelse(x > -500,
          -exp(dnorm(x, log = TRUE)
                - pnorm(x, log.p = TRUE))*x
           -exp(2*(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))),
          -1)
}
