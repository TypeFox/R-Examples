mkPYM <- function(PY,nyears) {
    PYM=.Call('mkPYM', PACKAGE = 'SEERaBomb', PY[,1], PY[,2], PY[,3],nyears)
#     print(nyears)
    dim(PYM)=c(126,nyears)
    colnames(PYM)=1973:(1973+nyears-1)
    rownames(PYM)=0.5:125.5
    PYM
}

