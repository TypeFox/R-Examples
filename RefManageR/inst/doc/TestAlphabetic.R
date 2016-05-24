## ----setup, include = TRUE, cache = FALSE--------------------------------
library(RefManageR)
bib <- ReadBib(system.file("Bib", "biblatexExamples.bib", 
                           package = "RefManageR"), check = FALSE)
BibOptions(check.entries = FALSE, style = "markdown", bib.style = "alphabetic", cite.style = 'alphabetic')

## ----fig.width=7, fig.height=6-------------------------------------------
plot(cars)

## ----fig.width=7, fig.height=6-------------------------------------------
plot(cars)

## ----results = "asis", echo = FALSE--------------------------------------
PrintBibliography(bib, .opts = list(check.entries = FALSE))

