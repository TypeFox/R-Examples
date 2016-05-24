library(catnet)

data("breast")

## first two columns and the first row are annotaions, skip them
psamples <- breast[-1, -c(1,2)]
## set the gene names
rownames(psamples) <- breast[-1, 1]
## convert the breast frame in gene-column format by transposing it
psamples <- as.data.frame(t(psamples))
## categorize the sample
psamples <- cnDiscretize(psamples, 2)

## select a subset of 100 genes
subset <- seq(2, 1201, 12)
psamples <- psamples[,subset]

eval <- cnSearchSA(data=psamples, perturbations=NULL, maxParentSet=1, maxComplexity=250, 
	parentsPool=NULL, fixedParents=NULL, 
	tempStart=100, tempCoolFact=0.75, tempCheckOrders=10, orderShuffles=4, stopDiff=0.000000001, maxIter=100, echo=TRUE)
eval

## select a network and plot it
cnet1 <- cnFind(eval, 200)
cnPlot(cnet1)
