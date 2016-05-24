## ----echo=FALSE----------------------------------------------------------
knitr::opts_chunk$set(dev="svg", fig.width=7, fig.height=7)

## ------------------------------------------------------------------------
library(BeviMed)
set.seed(0)

## ------------------------------------------------------------------------
y <- c(rep(TRUE, 20), rep(FALSE, 80))
Z <- c(rep(TRUE, 3), rep(FALSE, 20))

## ------------------------------------------------------------------------
G1 <- sapply(y, function(y_i) as.integer(runif(n=length(Z)) < 0.15))

probability_pathogenic(G=G1, y=y)

## ------------------------------------------------------------------------
G2 <- sapply(y, function(y_i) as.integer(runif(n=length(Z)) < 0.15 | 
	(if (y_i) 1:length(Z) == sample(which(Z), size=1) else rep(FALSE, length(Z)))))

probability_pathogenic(G=G2, y=y)

## ------------------------------------------------------------------------
output <- summary(bevimed(G=G2, y=y))

output

## ------------------------------------------------------------------------
print(output, print_Z=TRUE)

