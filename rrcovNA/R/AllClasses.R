## Define class unions for optional slots, e.g. for definition
##  of slots which will be computed on demand, like the
##  mahalanobis/robust distances
##setClassUnion("Uvector", c("vector", "NULL"))
##setClassUnion("Umatrix", c("matrix", "NULL"))
##setClassUnion("Ulist", c("list", "NULL"))

setClass("CovNA", representation("VIRTUAL"),
                    contains="Cov")

setClass("SummaryCovNA", representation(),
                    contains="SummaryCov")

setClass("CovNARobust", representation("VIRTUAL"),
                    contains=c("CovNA", "CovRobust"))

setClass("SummaryCovNARobust", representation(),
                    contains="SummaryCovRobust")

setClass("CovNAClassic", representation(),
                    contains=c("CovNA", "CovClassic"))

setClass("CovNAMcd", representation(),
                    contains=c("CovNARobust", "CovMcd"))

setClass("CovNAOgk", representation(),
                    contains=c("CovNARobust", "CovOgk"))

setClass("CovNASde", representation(),
                    contains=c("CovNARobust", "CovSde"))

setClass("CovNASest", representation(),
                    contains=c("CovNARobust", "CovSest"))

setClass("PcaNA", representation(Ximp = "Umatrix"), contains="Pca")

## setGeneric("impute", function(obj, ...) standardGeneric("impute", ...))
