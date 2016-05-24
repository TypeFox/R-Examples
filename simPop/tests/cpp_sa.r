rm(list=ls())
library(Rcpp)
library(simPop)
data(eusilcP)
data(eusilcS)

## synth pop:
pop <- eusilcP
colnames(pop)[3] <- "hhsize"

## donor data:
donors <- sampHH(pop, strata="region", hsize="hhsize")

## microdata (1) and donor (0)
pop <- rbind(pop, donors)
pop$weights <- rep(c(1,0), c(nrow(eusilcP), nrow(donors)))
pop$new.weights <- pop$weights

index <- pop$weights == 1
tab <- table(pop$region[index],pop$gender[index], pop$hhsize[index])

## create target marginals
totals <- tableWt(eusilcS[,c("db040", "rb090","hsize")], weights=eusilcS$rb050)
totals <- nrow(pop[index])/sum(totals) * totals
totals <- ceiling(totals)
totals <- as.data.table(totals)
setnames(totals, c("region","gender","hhsize","N"))

data <- data.table(pop, keep.rownames=TRUE)
totals <- totals
hid <- "hid"
parameter <- c("gender","hhsize")
split <- "region"
temp <- 10
eps.factor <-  0.02
maxiter <- 250
temp.cooldown <- 0.90
factor.cooldown <- 0.95
min.temp <- 10^-2
sample <- TRUE
parallel <- TRUE
verbose <- FALSE

sourceCpp("src/calibPop.cpp") 
source("R/calibPop.R")

data <- calibPop_cpp(data, totals, hid, parameter, split, temp = temp, eps.factor = eps.factor, 
                     maxiter=maxiter, temp.cooldown = temp.cooldown, factor.cooldown = factor.cooldown, min.temp = min.temp, verbose=verbose, parallel=parallel)

# compare with fixed margins
index <- data$weights == 1
tab_old <- as.data.table(table(data$region[index],data$gender[index], data$hhsize[index]))
setnames(tab_old, c("region",parameter,"before"))
setkeyv(tab_old,c("region",parameter))

index <- data$new.weights == 1
tab <- as.data.table(table(data$region[index],data$gender[index], data$hhsize[index]))
setnames(tab, c("region",parameter,"after"))
setkeyv(tab,c("region",parameter))

setnames(totals, c("region",parameter,"goal"))
setkeyv(totals,c("region",parameter))

res <- totals[tab_old]
res <- res[tab]
res
