library(catnet)

data("novartis")

## first two columns are annotaions, skip them
## take the first 100 genes only
psamples <- novartis[1:min(100,dim(novartis)[1]), 3:dim(novartis)[2]]
rownames(psamples) <- novartis[1:dim(psamples)[1],1]
## convert the novartis data frame in gene-column format by transposing it
psamples <- as.data.frame(t(psamples))
## categorize the sample
psamples <- cnDiscretize(psamples, 2)

eval <- cnSearchSA(data=psamples, perturbations=NULL, maxParentSet=2, maxComplexity=150, 
	parentsPool=NULL, fixedParents=NULL, 
	tempStart=1.0, tempCoolFact=0.75, tempCheckOrders=10, maxIter=100, echo=TRUE)
eval

## select a network and plot it
cnet1 <- cnFind(eval, 110)
cnPlot(cnet1)

