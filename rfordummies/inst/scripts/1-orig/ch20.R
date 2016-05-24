# Chapter 20 - Ten Tips on Working with Packages

## Poking Around the Nooks and Crannies of CRAN

options("repos" = c(CRAN = "http://cran.ma.imperial.ac.uk/"))

## Finding Interesting Packages

## Installing Packages

install.packages("fortunes")

## Loading Packages

library("fortunes")

## Reading the Package Manual and Vignette

library(help=fortunes)
vignette("fortunes")

## Updating Packages

update.packages()

## Unloading Packages

search()
detach(package:fortunes, unload=TRUE)

## Forging Ahead with R-Forge

install.packages("data.table", repos="http://R-Forge.R-project.org")

## Conducting Installations from BioConductor

source("http://bioconductor.org/biocLite.R")

## Reading the R Manual

