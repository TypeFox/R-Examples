## ---- echo=FALSE, results='asis'-----------------------------------------
data(dataKm, package = "RcmdrPlugin.KMggplot2")
exampleData <- dataKm[, c("time", "event", "trt", "sex", "marker")]
knitr::kable(head(exampleData, 10))

