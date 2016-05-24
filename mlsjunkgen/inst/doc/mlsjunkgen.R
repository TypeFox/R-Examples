## ----library, echo=FALSE-------------------------------------------------
library(mlsjunkgen)

## ----junkgen-------------------------------------------------------------
w <- 1
x <- 2
y <- 3
z <- 4
junkgen(w = w, x = x, y = y, z = z)

## ----mlsjunkgenv 1-------------------------------------------------------
mlsjunkgenv(n = 10, w = w, x = x, y = y, z = z, round = 8)

## ----mlsjunkgenv 2-------------------------------------------------------
mlsjunkgenv(n = 10, w = w, x = x, y = y, z = z)

## ----mlsjunkgend 1-------------------------------------------------------
mlsjunkgend(n = 10, w = w, x = x, y = y, z = z, round = 8)

## ----mlsjunkgend 2-------------------------------------------------------
mlsjunkgend(n = 10, w = w, x = x, y = y, z = z)

## ----mlsjunkgenm---------------------------------------------------------
mlsjunkgenm(nrow = 5, ncol = 5, w = w, x = x, y = y, z = z, round = 3)

