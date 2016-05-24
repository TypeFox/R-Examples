### R code from vignette source 'fortunes.Rnw'

###################################################
### code chunk number 1: fortunes.Rnw:19-25
###################################################
library("fortunes")
library("utils")
f <- read.fortunes(system.file("fortunes", "fortunes.csv", package = "fortunes"))
n <- nrow(f)
fortunes <- lapply(1:n, function(i) {toLatex(fortune(i, fortunes.data = f), number = TRUE, width = c(1, 0.85))})
invisible(lapply(fortunes, function(i) { print(i); cat("\\\\[0.2cm]\n\n") }))


