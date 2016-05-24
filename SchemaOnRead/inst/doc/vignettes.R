## ------------------------------------------------------------------------
library(SchemaOnRead)
results <- schemaOnRead("../inst/extdata")
print(results$dir1$Data.csv)

## ------------------------------------------------------------------------
library(SchemaOnRead)
xmlFile <- schemaOnRead("../inst/extdata/data.xml")
print(xmlFile)

## ------------------------------------------------------------------------
library(SchemaOnRead)
folder <- schemaOnRead("../inst/extdata", verbose = TRUE)

