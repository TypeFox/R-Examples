### R code from vignette source 'preface.Rnw'
### Encoding: ASCII

###################################################
### code chunk number 1: setup
###################################################
rm(list = ls())
if (!file.exists("tables")) dir.create("tables")
if (!file.exists("figures")) dir.create("figures")
set.seed(290875)
options(prompt = "R> ", continue = "+  ",
    width = 63, # digits = 4, 
    SweaveHooks = list(leftpar = function() 
        par(mai = par("mai") * c(1, 1.05, 1, 1))))
HSAURpkg <- require("HSAUR")
if (!HSAURpkg) stop("cannot load package ", sQuote("HSAUR"))
rm(HSAURpkg)
a <- Sys.setlocale("LC_ALL", "C")
book <- TRUE
refs <- cbind(c("AItR", "SI", "CI", "ANOVA", "MLR", "GLM", 
                "DE", "RP", "SA", "ALDI", "ALDII", "MA", "PCA", 
                "MDS", "CA"), 1:15)
ch <- function(x, book = TRUE) {
    ch <- refs[which(refs[,1] == x),]
    if (book) {
        return(paste("Chapter~\\\\ref{", ch[1], "}", sep = ""))
    } else {
        return(paste("Chapter~\\\\ref{", ch[2], "}", sep = ""))
    }
}


###################################################
### code chunk number 2: pre-HSAUR-install (eval = FALSE)
###################################################
## install.packages("HSAUR")


###################################################
### code chunk number 3: pre-HSAUR
###################################################
library("HSAUR")


###################################################
### code chunk number 4: pre-vignette (eval = FALSE)
###################################################
## vignette("Ch_introduction_to_R", package = "HSAUR")


###################################################
### code chunk number 5: pre-vignette-source (eval = FALSE)
###################################################
## edit(vignette("Ch_introduction_to_R", package = "HSAUR"))


###################################################
### code chunk number 6: pre-vignette-source (eval = FALSE)
###################################################
## vignette(package = "HSAUR")


