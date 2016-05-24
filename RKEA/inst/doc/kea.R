### R code from vignette source 'kea.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: kea.Rnw:21-24
###################################################
options(width = 75)
### for sampling
set.seed <- 1234


###################################################
### code chunk number 2: kea.Rnw:55-56
###################################################
library("RKEA")


###################################################
### code chunk number 3: kea.Rnw:65-80
###################################################
library("tm")
data("crude")

keywords <- list(c("Diamond", "crude oil", "price"),
                 c("OPEC", "oil", "price"),
                 c("Texaco", "oil", "price", "decrease"),
                 c("Marathon Petroleum", "crude", "decrease"),
                 c("Houston Oil", "revenues", "decrease"),
                 c("Kuwait", "OPEC", "quota"))

tmpdir <- tempfile()
dir.create(tmpdir)
model <- file.path(tmpdir, "crudeModel")

createModel(crude[1:6], keywords, model)


###################################################
### code chunk number 4: kea.Rnw:97-100
###################################################
extractKeywords(crude, model)

unlink(tmpdir, recursive = TRUE)


###################################################
### code chunk number 5: kea.Rnw:116-125 (eval = FALSE)
###################################################
## txts <- Sys.glob(file.path("fao780", "*.txt"))
## keys <- sub("txt$", "key", txts)
## txts <- lapply(txts, readLines)
## keys <- lapply(keys, readLines)
## build <- seq_len(100)
## xtrct <- seq(101, 105)
## model <- "fao780_model"
## createModel(txts[build], keys[build], model, "agrovoc", "skos")
## extractKeywords(txts[xtrct], model, "agrovoc", "skos")


