## ----c0, include = FALSE-------------------------------------------------
## rmarkdown::tufte_handout
library(knitr)
opts_chunk$set(comment = "", warning = FALSE, message = FALSE, fig.height = 4.94, 
               fig.width = 8, out.width = "690px", out.height = "426px")
oldOpt <- options(scipen = 10, width = 75)
library(dplyr)
library(fractional)
library(ggplot2)

## ------------------------------------------------------------------------
F <- c(1, 1, numeric(15))
for(i in 3:17) ## Fibonacci sequence by recurrence
  F[i] <- F[i-1] + F[i-2]  
F
vfractional((sqrt(5) + 1)/2, eps = 0, maxConv = 1:16)

## ----pdenom--------------------------------------------------------------
partial_denominators <- function(x, k = 10) {
  b <- rep(NA, k)
  r <- x
  for(i in 1:k) {
    b[i] <- floor(r)
    r <- r - b[i]
    if(isTRUE(all.equal(r, 0))) break
    r <- 1/r
  }
  structure(b, names = paste0("b", 1:k-1))
}

## ----pdenom_2, results = "asis"------------------------------------------
x <- c(pi = base::pi, e = exp(1), phi = (sqrt(5) + 1)/2, 
       structure(sqrt(1:9), names = paste0("sqrt(", 1:9, ")")))
tab <- x %>% sapply(partial_denominators) %>% t 
tab[is.na(tab)] <- ""
kable(tab, align = "r", caption = "Partial denominators")

## ----check, results = "hold"---------------------------------------------
pq <- .ratr(pi)
cat("Pn = ", pq[1], ", Qn = ", pq[2], ", n = ", pq[3], "\n", sep = "")
cat("pi = ", format(pi, digits = 15), 
    ", Pn/Qn = ", format(pq[1]/pq[2], digits = 15), 
    ", Error = ", pi - pq[1]/pq[2], "\n", sep = "")

## ----pi, echo = FALSE----------------------------------------------------
oldOpt <- options(scipen = 15, digits = 15)
pi_approx <- vfractional(base::pi, eps = 0, maxConv = 1:10)
tab <- within(data.frame(fraction = pi_approx, stringsAsFactors = FALSE), {
  value = numerical(fraction)
#  pi = base::pi
  error = base::pi - value
  n = seq_along(value) - 1
})[, c("n", "fraction", "value", "error")]
names(tab)[2] <- "Pn/Qn"
kable(tab, align = c("c", "c", "c", "r"))
options(oldOpt)

## ----sqrt5---------------------------------------------------------------
(s5 <- vfractional(sqrt(5), eps = 0, maxConv = 1:7))
d5 <- denominators(s5)
e5 <- abs(sqrt(5) - numerical(s5))

## ----sqrt5_2-------------------------------------------------------------
d <- seq(max(d5))
n <- round(sqrt(5) * d)

## ------------------------------------------------------------------------
gcd <- mapply(FUN = function(a, b) if(b == 0) a else Recall(b, a %% b),
              n, d)
nd <- cbind(n, d)/gcd
nd <- nd[!duplicated(nd), ]
e <- abs(sqrt(5) - nd[, 1]/nd[, 2])

## ----sqrt5_3,echo=FALSE--------------------------------------------------
dat <- data.frame(Denominator = nd[,2], Error = e)
dat5 <- data.frame(Denominator = d5, Error = e5)
ggplot(dat) + aes(x = Denominator, y = Error) + 
  geom_point(colour = "steel blue", size = 0.5) +
  scale_x_log10(breaks = 5^(0:5)) + scale_y_log10() +
  geom_point(data = dat5, mapping = aes(x = Denominator, y = Error), colour = "red",
             shape = 1, size = 3) +
  geom_step(data = dat5, mapping = aes(x = Denominator, y = Error),
            colour = "red", size = 0.5) +
  xlab(expression(paste("Denominator, ", italic(d)))) +
  ylab("Absolute error") +
  ggtitle(expression(paste("Errors in Rational Approximations, ",
                             italic(n)/italic(d), ", to ", sqrt(5))))

