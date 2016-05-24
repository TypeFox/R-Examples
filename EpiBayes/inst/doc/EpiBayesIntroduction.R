## ----setup, include=FALSE------------------------------------------------
library(knitr)
#opts_chunk$set(cache=TRUE, autodep=TRUE, fig.width=5, fig.height=5)

## ----package_loading, hide=TRUE, error=FALSE, warning=FALSE, message=FALSE----
library(epiR)  # For the BetaBuster function
library(compiler)  # To compile the larger functions for computational speed
library(coda)  # For processing Bayesian model output
library(shape)  # For nice colorbar legends
library(scales)  # For transparent colors

## ----model_source, eval = FALSE------------------------------------------
#  install.packages("EpiBayes_0.0.1", type = "source", repos = NULL)  # Make sure the
#  	version is correct and the working directory is pointed to where the .tar.gz
#  	file is stored
#  library(EpiBayes)  # Load the package

## ----eval=FALSE----------------------------------------------------------
#  name_of_your_model$taumat[1, 1, ]

## ----eval=FALSE----------------------------------------------------------
#  name_of_your_model$taumat[2, 1, ]

## ----eval=FALSE----------------------------------------------------------
#  hist(name_of_your_model$taumat[1, 1, ], col = "cyan");box("plot")

## ----eval=FALSE----------------------------------------------------------
#  plot(name_of_your_model$taumat[1, 1, ], type = "l")

## ----eval=FALSE----------------------------------------------------------
#  plot(name_of_your_model$taumat[1, 1, -c(1:1000)], type = "l")

## ----eval=FALSE----------------------------------------------------------
#  name_of_your_model$pimat[10, 1, ]

