## ----setup, include = FALSE----------------------------------------------
library(qlcMatrix)

## ----ttIntro-------------------------------------------------------------
data <- c("a", "b", "a", "c", "b", "c", "A")
ttMatrix(data)

## ----ttSimplify----------------------------------------------------------
ttMatrix(data, simplify = TRUE, collation.locale = "en_US.UTF-8")

## ----pwIntro-------------------------------------------------------------
strings <- c("this", "is", "that")
(PW <- pwMatrix(strings, sep = "", simplify = TRUE, gap.length = 1, gap.symbol = "_") )

## ----pwUsage-------------------------------------------------------------
TT <- ttMatrix(rownames(PW), simplify = TRUE)
printSpMatrix(TT, col.names = TRUE)
(TT*1) %*% (PW*1)

## ----adjacency-----------------------------------------------------------
S <- bandSparse( n = ncol(TT), k = 1)
(TT*1) %*% (S*1) %*% t(TT*1)

