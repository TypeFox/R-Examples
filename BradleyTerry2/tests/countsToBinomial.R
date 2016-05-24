library(BradleyTerry2)
data(citations, package = "BradleyTerry2")

## Convert frequencies to success/failure data
results <- countsToBinomial(citations)
results
