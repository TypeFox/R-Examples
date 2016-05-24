# ExpFeedClass
setClass("dataSet", representation(data = "matrix", binary = "logical", minimum = "numeric", maximum = "numeric", halfStar = "logical"))

# SVDclass####
setClass("SVDclass", representation(alg = "character", data = "dataSet", factors = "list", parameters = "list"))

# IBclass####
setClass("IBclass", representation(alg = "character", data = "dataSet", sim = "matrix", sim_index_kNN = "matrix", neigh = "numeric"))

# wALSclass####
setClass("wALSclass", representation(alg = "character", data = "dataSet", factors = "list", weightScheme = "matrix", parameters = "list"))

# BPRclass####
setClass("BPRclass", representation(alg = "character", data = "dataSet", factors = "list", parameters = "list"))

# PPLclass####
setClass("PPLclass", representation(alg = "character", data = "dataSet", indices = "numeric"))

# algAverageClass####
setClass("algAverageClass", representation(alg = "character", data = "dataSet", average = "matrix"))

# evalModel####
setClass("evalModel", representation(data = "dataSet", folds = "numeric", fold_indices = "list", fold_indices_x_user = "list"))

setClass("recResultsClass", representation(indices = "list", recommended = "list"))


# #evalResults#### setClass('evalResults', representation( )) 
