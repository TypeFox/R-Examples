## ------------------------------------------------------------------------
library(decode)

## ------------------------------------------------------------------------
path = system.file('extdata', package='decode')
geneSetInputFile = file.path(path, "geneSet.txt")
geneExpressionFile = file.path(path, "Expression_data_50genes.txt")
runDecode(geneSetInputFile, geneExpressionFile)


