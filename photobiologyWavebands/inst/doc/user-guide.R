## ----setup, include=FALSE, cache=FALSE---------------
library(knitr)
opts_chunk$set(fig.path='figure/pos-', fig.align='center', fig.show='hold', size="footnotesize",            fig.width=8, fig.height=5.5, out.width='.95\\textwidth')
options(replace.assign = TRUE, width = 55,
        warnPartialMatchAttr = FALSE,
        warnPartialMatchDollar = FALSE,
        warnPartialMatchArgs = FALSE)

## ----example-0-hiden, include=FALSE------------------
library(photobiology)
library(photobiologyWavebands)

## ----own-set-up, echo=FALSE, include=FALSE-----------
my_version <- packageVersion("photobiologyWavebands")

## ----example-1---------------------------------------
e_irrad(sun.spct, UV()) # W m-2
q_irrad(sun.spct, UV()) * 1e6 # umol s-1 m-2

## ----------------------------------------------------
e_irrad(sun.spct, list(Blue(), VIS()))
e_irrad(sun.spct, list(B = Blue(), VIS()))

## ----------------------------------------------------
e_irrad(sun.spct, VIS_bands())

## ----------------------------------------------------
q_ratio(sun.spct, Blue(), VIS())

## ----------------------------------------------------
e_irrad(sun.spct, CIE())

## ----eval=TRUE---------------------------------------
# at 1 nm intervals
wavelengths1 <- 285:400
action.spectrum1 <- CIE_e_fun(wavelengths1)

## ----------------------------------------------------
sun.spct * CIE()

## ----------------------------------------------------
e_response(sun.spct * CIE1924_lef.spct) * photopic_sensitivity

## ----------------------------------------------------
e_response(sun.spct * CIE2008_lef2deg.spct) * photopic_sensitivity

## ----------------------------------------------------
e_response(sun.spct * CIE2008_lef2deg.spct) * photopic_sensitivity *
                       interpolate_spct(CIE2008_lef2deg.spct, 555)$s.e.response

## ----------------------------------------------------
e_response(sun.spct * 1e-6 * CIE1951_scotopic_lef.spct) * scotopic_sensitivity

