## ----ini, echo=FALSE, results='hide'-------------------------------------
library(knitr)

## ----echo=FALSE----------------------------------------------------------
opts_chunk$set(
    echo = FALSE,
    comment = NA,
    quiet = TRUE,
    progress = FALSE,
    tidy = FALSE,
    cache = FALSE,
    message = FALSE,
    error = TRUE,
    warning = TRUE
)

## ----echo=FALSE----------------------------------------------------------
rversion <- sub(
    pattern = "R *\\(>= *([^)]*)\\).*", 
    replacement = "\\1", 
    x = packageDescription("BEQI2", fields = "Depends")
)

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  library(BEQI2)

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2dir()

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2dir("myBEQI2_analysis_dir")

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  BEQI2()

## ----eval=FALSE, echo=TRUE-----------------------------------------------
#  beqi2()

