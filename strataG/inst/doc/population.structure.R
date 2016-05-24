## ----echo = FALSE, message = FALSE---------------------------------------
options(digits = 2)
library(strataG)

## ------------------------------------------------------------------------
data(msats.g)
msats <- stratify(msats.g, "fine")
msats <- msats[, locNames(msats)[1:4], ]

## ------------------------------------------------------------------------
statFst(msats)

statGst(msats, nrep = 10, keep.null = TRUE)

## ------------------------------------------------------------------------
ovl <- overallTest(msats, stats = c("Fst", "Chi2"), nrep = 1000)

## ------------------------------------------------------------------------
pws <- pairwiseTest(msats, stats = c(statFstPrime, statGst), nrep = 1000)

## ------------------------------------------------------------------------
pws

## ------------------------------------------------------------------------
popStruct <- popStructTest(msats, stats = c(statFst, statFstPrime), nrep = 1000, quietly = TRUE)
popStruct

