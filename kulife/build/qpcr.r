# Hent read.xls

library(gdata)

indata <- read.xls("qpcr.xls", header=TRUE)


cycle <- rep(1:45, 162/3)

line  <- rep(c("wt", "rnt", "empty"), times=45*c(20, 20, 14))


flour <- indata[, seq(3, 162, 3)]

qpcr <- data.frame(flour= as.vector(as.matrix(flour[,c(1:7,seq(21,33,2))])),
                   lines=factor(rep(c("wt", "rnt"), times=45*c(7, 7))),
                   cycle=rep(1:45, 14),
                   experiment=rep(1:14,times=rep(45,14)))

