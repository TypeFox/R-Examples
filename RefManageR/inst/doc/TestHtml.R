## ----setup, include = TRUE, cache = FALSE--------------------------------
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib",
                           package = "RefManageR"), check = FALSE)
bib2 <- ReadBib(system.file("Bib", "RJC.bib", package = "RefManageR"))[[seq_len(20)]]
BibOptions(check.entries = FALSE, style = "html")

## ----fig.width=7, fig.height=6-------------------------------------------
plot(cars)

## ----fig.width=7, fig.height=6-------------------------------------------
plot(cars)

## ----results = "asis", echo = FALSE--------------------------------------
PrintBibliography(bib, .opts = list(check.entries = FALSE, bib.style = "draft"))

## ----results = "asis", echo = FALSE--------------------------------------
PrintBibliography(bib2)

