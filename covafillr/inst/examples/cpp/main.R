
library(covafillr)
require(TMB)

## From TMB
## https://github.com/kaskr/adcomp/blob/eb9268a2865905023e14bcadf774aa4fae7508b3/TMB/R/TMB.R


compiler <- c("g++","clang++")

ppflags <- paste(cxxFlags(),
                 paste0("-I",system.file("include",package="RcppEigen"))
                 )


d <- list()

for(i in 1:length(compiler)){
a <- system2(compiler[i],
        args = c(cxxFlags(),
                 paste0("-I",system.file("include",package="RcppEigen")),
                 'test.cpp'),
        stdout=TRUE,stderr=TRUE
        )
d[[i]] <- c(system2('./a.out',stdout=TRUE,stderr=TRUE))
file.remove('a.out')
}
