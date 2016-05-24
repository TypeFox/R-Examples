library(matconv)
dataConvs <- list()
dataConvs[[1]] <- makeDataMap(matClass = "string", rClass = "vector")
dataConvs[[2]] <- makeDataMap(matClass = "matrix", rClass = "matrix")
dataConvs[[2]] <- makeDataMap("{", "}", "matrix")

dataConvs[[3]] <- makeSliceMap("{", "}", "list")
dataConvs[[4]] <- makeSliceMap(matClass = "structure", rClass = "list")
