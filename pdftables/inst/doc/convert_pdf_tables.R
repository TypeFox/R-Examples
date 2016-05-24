## ----options, echo=FALSE-------------------------------------------------
knitr::opts_chunk$set(cache = FALSE, warning = FALSE, error = FALSE, 
                      fig.path = "", eval=FALSE)

## ----install-------------------------------------------------------------
#  # From CRAN
#  install.packages("pdftables")
#  
#  # From Github
#  library(devtools)
#  install_github("expersso/pdftables")
#  
#  library(pdftables)

## ------------------------------------------------------------------------
#  write.csv(head(iris, 20), file = "test.csv", row.names = FALSE)
#  
#  # Open test.csv and print as PDF to "test.pdf"
#  
#  convert_pdf("test.pdf", "test2.csv")
#  # Converted test.pdf to test2.csv

## ------------------------------------------------------------------------
#  get_remaining()

