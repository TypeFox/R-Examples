## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(comment = "#>")

## ---- eval=TRUE, echo=FALSE----------------------------------------------
library(easypower)

## ---- eval = TRUE--------------------------------------------------------
# Define main effects
main.eff1 <- list(name = "A", levels = 2, eta.sq = 0.123)
main.eff2 <- list(name = "C", levels = 4, eta.sq = 0.215)

## ---- eval=FALSE---------------------------------------------------------
#  # Example of using the default eta.sq setting
#  main.eff <- list(name = "A", levels = 3)

## ---- eval=TRUE----------------------------------------------------------
# Change interaction effect size
int.eff1 <- list(name = "A*C", eta.sq = 0.079)

## ---- eval = TRUE, collapse=TRUE-----------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, int1 = int.eff1)

## ---- eval = TRUE, collapse=TRUE-----------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, interaction.eta2 = 0.079)

## ---- eval=TRUE, collapse=TRUE-------------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, int1 = int.eff1, result = "highest")

## ---- eval=TRUE, collapse=TRUE-------------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, int1 = int.eff1, result = "select")

## ---- eval = TRUE--------------------------------------------------------
# Define main effects
main.eff1 <- list(name = "Sex", levels = 2, eta.sq = 0.0099)
main.eff2 <- list(name = "Age", levels = 3, eta.sq = 0.0588)
main.eff3 <- list(name = "Conditions", levels = 4, eta.sq = 0.1506)

## ---- eval=FALSE---------------------------------------------------------
#  # Example of using the default eta.sq setting
#  main.eff <- list(name = "A", levels = 3)

## ---- eval=TRUE----------------------------------------------------------
# Changing the effect sizes of specific interactions
int.eff1 <- list(name = "Age*Conditions", eta.sq = "med")
int.eff2 <- list(name = "Sex*Conditions", eta.sq = "med")

## ---- eval=TRUE, collapse=TRUE-------------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, interaction.eta2 = 0.0588)

## ---- eval=TRUE, collapse=TRUE-------------------------------------------
n.multiway(iv1 = main.eff1, iv2 = main.eff2, iv3 = main.eff3, int1 = int.eff1, int2 = int.eff2)

