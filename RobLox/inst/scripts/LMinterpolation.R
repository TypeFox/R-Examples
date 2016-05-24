###############################################################################
## Interpolated functions to speed up computation of Lagrange Multipliers
## 
## regarding a change to .C calls and approxfun in R 2.16.0, we need to make
## distinction between version before 2.16.0 and afterwards
###############################################################################

library(RobLox)
radius <- c(1e-8, 5e-8, 1e-7, 5e-7, 1e-6, 5e-6, 1e-5, 5e-5, seq(1e-4, 0.01, by = 0.001),
            seq(0.02, 5, by = 0.01), seq(5.05, 10, by = 0.05))
location <- sapply(radius, rlOptIC, computeIC = FALSE)
scale <- sapply(radius, rsOptIC, computeIC = FALSE)

fun <- function(radius){
  print(radius)
  rlsOptIC.AL(radius, computeIC = FALSE)
}
locationScale <- sapply(radius, fun)
#locationScale <- sapply(radius, rlsOptIC.AL, computeIC = FALSE)

## location
.A.loc <- unlist(location[1,])
.b.loc <- unlist(location[3,])

## scale
.A.sc <- unlist(scale[1,])
.a.sc <- unlist(scale[2,])
.b.sc <- unlist(scale[3,])

## location and scale
n <- length(radius)
.A1.locsc <- unlist(locationScale[1,])[seq(1, 4*n-3, by = 4)]
.A2.locsc <- unlist(locationScale[1,])[seq(4, 4*n, by = 4)]
.a.locsc <- unlist(locationScale[2,])[seq(2, 2*n, by = 2)]
.b.locsc <- unlist(locationScale[3,])

.radius.gitter <- radius

## Saving the results in sysdata.rda
#load("sysdata.rda")
save(.radius.gitter, .finiteSampleRadius.loc, .finiteSampleRadius.locsc, .finiteSampleRadius.sc,
     .A.loc, .b.loc, .A.sc, .a.sc, .b.sc, .A1.locsc, .A2.locsc, .a.locsc, .b.locsc,
     file = "sysdata.rda")
