
## ----pkgmaker_preamble, message = FALSE, echo=FALSE, results='asis'------
library(pkgmaker)
pkgmaker::latex_preamble()


## ----hook_try------------------------------------------------------------
library(knitr)
knit_hooks$set(try = pkgmaker::hook_try)


## ----without_try---------------------------------------------------------
try( stop('this error will not appear in the document'))


## ----with_try, try = NA--------------------------------------------------
txt <- 'this error will be shown'
try( stop(txt) )


## ----with_try_highlight, try = TRUE--------------------------------------
txt <- 'this error will be shown'
try( stop(txt) )


## ----sessionInfo, echo=FALSE, results='asis'-----------------------------
toLatex(sessionInfo())


