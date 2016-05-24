## enforce par(ask=FALSE)
options(device.ask.default=FALSE)

## make some size reductions
tools::compactPDF(paths = "E:/R-Pkg-Produce/robustfa_1.0-5/vignettes/FaModel.pdf", gs_quality = "ebook")
tools::compactPDF(paths = "E:/R-Pkg-Produce/robustfa_1.0-5/vignettes/Structure.pdf", gs_quality = "ebook")

## create a LaTeX file "robustfa.tex"
Sweave("robustfa.Rnw")

## This can be compiled to PDF by
## tools::texi2pdf("robustfa.tex")
## or outside R by
## Rcmd texify --pdf robustfa.tex
## if MiKTeX is available.
tools::texi2pdf("robustfa.tex")

## make some size reductions
tools::compactPDF(paths = "E:/R-Pkg-Produce/robustfa_1.0-5/inst/doc/robustfa.pdf", 
qpdf = Sys.getenv("R_QPDF"), gs_cmd = Sys.getenv("R_GSCMD"), gs_quality = "ebook")


Sys.getenv("R_QPDF")
Sys.getenv("qpdf")
Sys.getenv("Path")
Sys.which(Sys.getenv("R_QPDF", "qpdf"))
Sys.getenv("R_GSCMD", "")
Sys.getenv("GS_QUALITY", "none")


## create an R source file from the code chunks and obtain "robustfa.R"
Stangle("robustfa.Rnw")

## which can be sourced, e.g.
source("robustfa.R")