library(testthat)
for(i in dir("../../R/")){
    source(paste("../../R/", i, sep = ''))
}
