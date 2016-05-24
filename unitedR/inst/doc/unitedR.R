## ----eval=FALSE----------------------------------------------------------
#  install.packages("unitedR")

## ----message=F-----------------------------------------------------------
library(unitedR)

## ------------------------------------------------------------------------
(home <- formation(10, NA, c(7,5,3), c(8,8), c(10,10,8,5,0)))
(away1 <- formation(5, 8, c(8,8,0,0), c(10,10), c(10,10,10),
 hardness = c(0,0,0,0,1)))
(away2 <- formation(10, 8, c(8,10), c(10,10), c(10,10,10,5,0),
 hardness = c(0,0,0,0,1), homeAdv = c(0,0,2,0,0)))
# unitedSim and unitedSimOne are similar in this particular case
unitedSim(home, away1)
unitedSim(home, away1, away2)

## ------------------------------------------------------------------------
set.seed(123)
(home <- formation(10, NA, c(7,5,3), c(8,8), c(10,10,8,5,0),
                   hardness = c(0,0,4,2,1)))
(away1 <- formation(5, 8, c(8,8,0,0), c(10,10), c(10,10,10),
 hardness = c(0,0,0,0,1)))
(away2 <- formation(10, 8, c(8,10), c(10,10), c(10,10,10,5,0),
 hardness = c(0,0,0,0,8), homeAdv = c(0,0,2,0,0)))
# unitedSim and unitedSimOne are similar in this particular case
unitedSim(home, away1, r = 100)
unitedSim(home, away1, away2, r = 100)

## ------------------------------------------------------------------------
(home <- formation(10, NA, 14, 14, 42))
(away1 <- formation(5, 8, 10, 10, 30))
(away2 <- formation(10, 8, 16, 16, 30, homeAdv = c(0,0,2,0,0)))
# unitedSim and unitedSimOne are similar in this particular case
unitedSim(home, away1)
unitedSim(home, away1, away2)

