## ----setup, include = FALSE----------------------------------------------
library(qlcData)

## ---- eval = FALSE-------------------------------------------------------
#  # install devtools from CRAN
#  install.packages("devtools")
#  # install qlcData from github using devtools
#  devtools::install_github("cysouw/qlcData")
#  # load qlcTokenize package
#  library(qlcData)
#  # access help files of the package
#  help(qlcData)

## ------------------------------------------------------------------------
test <- "hállo hállо"

## ---- eval = FALSE-------------------------------------------------------
#  write.profile(test)

## ----echo=FALSE, results='asis'------------------------------------------
# some example string
knitr::kable(write.profile(test))

## ------------------------------------------------------------------------
# the differenec between various "o" characters is mostly invisible on screen
"o" == "o"  # these are the same "o" characters, so this statement in true
"o" == "о"  # this is one latin and and cyrillic "o" character, so this statement is false

## ------------------------------------------------------------------------
test <- c("this thing", "is", "a", "vector", "with", "many", "strings")

## ---- eval = FALSE-------------------------------------------------------
#  write.profile(test)

## ----echo=FALSE, results='asis'------------------------------------------
# some example string
knitr::kable(write.profile(test))

## ------------------------------------------------------------------------
tokenize(test)

## ---- eval = FALSE-------------------------------------------------------
#  dir.create("~/Desktop/tokenize")
#  setwd("~/Desktop/tokenize")
#  tokenize(test, file.out = "test_profile.txt")

## ---- echo = FALSE, results='asis'---------------------------------------
test_profile.txt <- as.data.frame(rbind(as.matrix(tokenize(test)$profile),c(" ", "th"),c(" ","ng")))
knitr::kable(test_profile.txt)

## ---- eval = FALSE-------------------------------------------------------
#  tokenize(test, profile = "test_profile.txt")
#  
#  # with overwriting of the existing profile:
#  # tokenize(test, profile = "test_profile.txt", file.out = "test_profile.txt")
#  
#  # note that you can abbreviate this in R:
#  # tokenize_old(test, p = "test_profile.txt", f = "test_profile.txt")

## ---- echo = FALSE-------------------------------------------------------
tokenize(test, profile = test_profile.txt)

## ---- eval = FALSE-------------------------------------------------------
#  tokenize(c("think", "thin", "both"), profile = "test_profile.txt")

## ---- echo = FALSE-------------------------------------------------------
tokenize(c("think", "thin", "both"), profile = test_profile.txt)

## ---- echo = FALSE, results='asis'---------------------------------------
Grapheme <- c("c", "c", "n", "s", "a", "i")
IPA <- c("k", "tʃ", "n", "s", "a", "i")
Right <- c("", "[ie]", "", "", "", "")
italian <- cbind(Grapheme, Right, IPA)
knitr::kable(italian)

## ------------------------------------------------------------------------
tokenize(c("casa", "cina"), profile = italian, transliterate = "IPA", regex = TRUE)$strings

## ---- echo = FALSE, results='asis'---------------------------------------
Grapheme <- c("c", "c", "n", "s", "a", "i", "e")
IPA <- c("k", "tʃ", "n", "s", "a", "i", "e")
Right <- c("", "frontV", "", "", "", "","")
Class <- c("","","","","","frontV","frontV")
italian <- cbind(Grapheme, Right, Class, IPA)
knitr::kable(italian)

## ------------------------------------------------------------------------
tokenize(c("casa", "cina"), profile = italian, transliterate = "IPA", regex = TRUE)$strings

