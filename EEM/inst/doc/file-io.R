## ----setup, echo=FALSE, results="hide", warning=FALSE--------------------
suppressPackageStartupMessages({
    library(EEM)
})

## ----setup2, echo=FALSE, warning=FALSE-----------------------------------
library(knitr)
library(EEM)
opts_chunk$set(fig.path = "figure_file-io/")

## ----readEEM, eval=FALSE-------------------------------------------------
#  library(EEM)
#  
#  # importing files
#  data <-readEEM("sample1.txt") # read in a file
#  data <-readEEM(c("sample1.txt", "sample2.txt")) # read in two files
#  
#  # importing folders
#  data <-readEEM("/data") # read in all files in data folder
#  data <-readEEM(c("/data", "/data2")) # read in all files in two folders
#  data <-readEEM("C:\\data") # full path. Note that the slash is doubled.
#  data <- readEEM("C:/data") # read in all files in data folder. Aside from double slashes,
#                             # a reverted slash can also be used.

## ----import, eval=FALSE--------------------------------------------------
#  datamatrix <- read.csv("datamatrix.csv", row.names = 1) # for csv
#  datamatrix <- read.csv("datamatrix.txt", row.names = 1) # for txt
#  datamatrix[1:5,1:5] # check the first 5 rows and columns

## ----import_excel, eval=FALSE--------------------------------------------
#  library(readxl)
#  datamatrix <- read_excel("datamatrix.csv")
#  datamatrix[1:5,1:5] # check the first 5 rows and columns

## ----import_excel2, eval=FALSE-------------------------------------------
#  rownames(datamatrix) <- datamatrix[,1] # assign the first column to the matrix row names
#  datamatrix <- datamatrix[,-1] # delete the first column
#  datamatrix[1:5,1:5] # check the first 5 rows and columns

## ----fold, eval=FALSE----------------------------------------------------
#  library(EEM)
#  datamatrix_folded <- fold(as.matrix(datamatrix))
#  class(datamatrix_folded)
#  # "EEM"

