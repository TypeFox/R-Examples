## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, prompt=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  externaldata <- read.csv("externaldata.csv", header=TRUE)
#  externaldata$match <- fullmatch(..., data=externaldata)
#  write.csv(externaldata, file="externaldata.matched.csv")

## ----eval=FALSE----------------------------------------------------------
#  sasdata <- read.csv("C:/Users/myuser/Desktop/sasout.csv", header=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  sasdata <- read.csv("C:\\Users\\myuser\\Desktop\\sasout.csv", header=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(optmatch)
#  f <- fullmatch(gender ~ age + ppty, data=sasdata)
#  sasdata$match <- f

## ----eval=FALSE----------------------------------------------------------
#  write.csv(sasdata, "C:/Users/myuser/Desktop/rout.sas.csv", row.names=FALSE)

## ----eval=FALSE----------------------------------------------------------
#  library(foreign)
#  write.foreign(sasdata, 'C:/Users/myuser/Desktop/rout.sas.txt',
#                'C:/Users/myuser/Desktop/rout.code.sas', package = 'SAS')

## ----eval=FALSE----------------------------------------------------------
#  statadata <- read.csv("C:/Users/myuser/Desktop/stataout.csv", header=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  statadata <- read.csv("C:\\Users\\myuser\\Desktop\\stataout.csv", header=TRUE)

## ----eval=FALSE----------------------------------------------------------
#  library(optmatch)
#  f <- fullmatch(foreign ~ price + mpg + ppty, data=statadata, max.controls=3)
#  statadata$match <- f

## ----eval=FALSE----------------------------------------------------------
#  write.dta(statadata, "C:/Users/myuser/Desktop/rout.stata.dta")

