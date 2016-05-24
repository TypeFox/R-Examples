## ----c0, include = FALSE---------------------------------------------------------------------
## rmarkdown::tufte_handout
library(knitr)
opts_chunk$set(comment = "", warning = FALSE, message = FALSE, fig.height = 4.94, 
               fig.width = 8, out.width = "690px", out.height = "426px")
oldOpt <- options(scipen = 10, width = 95)

## ----c1, results = "asis", echo = FALSE------------------------------------------------------
Hilbert <- function(n) {
  structure(outer(1:n, 1:n, function(r, c) 1/(r+c-1)), 
            dimnames = list(paste0("R",1:n), paste0("C", 1:n)))
}
H5 <- Hilbert(5)
kable(H5, align = "c")

## ----c2, results = "asis", echo = FALSE------------------------------------------------------
library(dplyr)
library(fractional)
H5 %>% fractional %>% as.character %>% kable(align = "r")

## --------------------------------------------------------------------------------------------
library(dplyr)
(x <- matrix(1:9/12, 3, 3) %>% fractional)
(xr <- 1/x)
x + 10
(x2 <- x^2)
(sx2 <- sqrt(x2))  ## demoted back to numeric by sqrt()
fractional(sx2)    ## reinstated

solve(xr) %>% fractional  ## xr is non-singular with a simple inverse

numerators(x2)     ## numerators and denominators are available
denominators(x2)

## ---- include = FALSE---------------------------------------------------------
oldWidth <- options(width = 80)

## -----------------------------------------------------------------------------
F <- c(1,1,numeric(15))
for(i in 3:17) ## Fibonacci sequence by recurrence
  F[i] <- F[i-1] + F[i-2]  
F
(phi <- (sqrt(5) + 1)/2)

## -----------------------------------------------------------------------------
vfractional(phi, eps = 0, maxConv = 1:16)

## ---- include = FALSE------------------------------------------------------------------------
options(oldWidth)

## --------------------------------------------------------------------------------------------
library(ggplot2)
N <- 500000
set.seed(3210)

rat(rnorm(N), eps = 1.0e-07)[, "n"] %>% 
  table(x = .)                      %>% 
  as.data.frame(responseName = "f") %>% 
  ggplot(.) + aes(x = x, y = f/N) +
  geom_bar(stat = "identity", fill = "steel blue") +
  ylab("Relative frequency") + xlab("Convergents")

