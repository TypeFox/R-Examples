## ----setup, include=FALSE, cache=FALSE-----------------------------------
library("knitr")
opts_chunk$set(tidy.opts=list(width.cutoff=45), tidy=FALSE, fig.align="center",
               fig.height=4.25, comment=NA, prompt=TRUE)

## ----env, echo=FALSE-----------------------------------------------------
suppressPackageStartupMessages(library("MALDIquant"))
suppressPackageStartupMessages(library("MALDIquantForeign"))

## ----fileformats---------------------------------------------------------
supportedFileFormats()

## ----mqsetup, eval=FALSE-------------------------------------------------
#  install.packages(c("MALDIquant", "MALDIquantForeign"))

## ----mqlibrary, eval=FALSE-----------------------------------------------
#  library("MALDIquant")
#  library("MALDIquantForeign")

## ----import--------------------------------------------------------------
## get the example directory
exampleDirectory <- system.file("exampledata",
                                package="MALDIquantForeign")

spectra <- import(file.path(exampleDirectory,
                            "brukerflex"),
                  verbose=FALSE)
spectra[[1]]

## ----importbrukerflex----------------------------------------------------
spectra <- importBrukerFlex(file.path(exampleDirectory,
                                      "brukerflex"),
                            verbose=FALSE)
spectra[[1]]

## ----importcsvcompressed-------------------------------------------------
spectra <- importCsv(file.path(exampleDirectory, "compressed",
                               "csv.tar.gz"), verbose=FALSE)
spectra[[1]]

spectra <- importCsv(file.path(exampleDirectory, "compressed",
                               "csv.zip"), verbose=FALSE)
spectra[[1]]

## ----importremote, eval=FALSE--------------------------------------------
#  spectra <- import(paste0("http://www.meb.ki.se/",
#                           "~yudpaw/papers/spikein_xml.zip"),
#                    centroided=FALSE, verbose=TRUE)

## ----centroided----------------------------------------------------------
peaks <- import(file.path(exampleDirectory, "ascii.txt"),
                centroided=TRUE, verbose=FALSE)
peaks

## ----masspectrum---------------------------------------------------------
spectra <- list(
  createMassSpectrum(mass=1:5, intensity=1:5),
  createMassSpectrum(mass=1:5, intensity=6:10))

## ----export1-------------------------------------------------------------
export(spectra[[1]], file="spectrum1.csv")
import("spectrum1.csv")

## ----exportpath, eval=TRUE-----------------------------------------------
export(spectra, type="csv", path="spectra", force=TRUE)
list.files("spectra")

## ----mqvignette, eval=FALSE----------------------------------------------
#  vignette(topic="MALDIquant", package="MALDIquant")

## ----sessioninfo, echo=FALSE, results="asis"-----------------------------
toLatex(sessionInfo(), locale=FALSE)

