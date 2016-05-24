## enforce par(ask=FALSE)
options(device.ask.default=FALSE)

## create a LaTeX file "robustfa.tex"
Sweave("OneTwoSamples.Rnw")

## This can be compiled to PDF by
## tools::texi2pdf("OneTwoSamples.tex")
## or outside R by
## Rcmd texify --pdf OneTwoSamples.tex
## if MiKTeX is available.
tools::texi2pdf("OneTwoSamples.tex")

## create an R source file from the code chunks and obtain "robustfa.R"
Stangle("OneTwoSamples.Rnw")

## which can be sourced, e.g.
source("OneTwoSamples.R")