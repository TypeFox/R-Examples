setClass("Outlier", representation(call = "language",
                              counts = "numeric",
                              grp = "factor",
                              wt = "Uvector",
                              flag = "Uvector",
                              method = "character",
                              singularity = "Ulist",
                              "VIRTUAL")
)

setClass("OutlierMahdist", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierPCOut", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierPCDist", representation(covobj="Ulist",
                                    k="numeric"),
                                    contains="Outlier"
)

setClass("OutlierSign1", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("OutlierSign2", representation(covobj="Ulist"),
                                    contains="Outlier"
)

setClass("SPcaGrid", representation(),
                    contains="PcaGrid")

###################### SIMCA ####################################
setClass("Simca", representation(call = "language",
                               prior = "vector",
                               counts = "vector",
                               pcaobj="Ulist",
                               k = "Uvector",
                               flag = "Uvector",
                               X = "Umatrix",
                               grp = "factor",
                               "VIRTUAL"))

setClass("CSimca", contains="Simca")
setClass("RSimca",  contains="Simca")
setClass("PredictSimca", representation(classification = "factor",
                                      odsc = "matrix",
                                      sdsc = "matrix",
                                      ct="Utable"))
setClass("SummarySimca", representation(simcaobj = "Simca"))
