
## This chisq test is used in the Leong et al paper

library(RcppZiggurat)

bins <- 200
chires <- RcppZiggurat:::chisqTest(draws=1e8,          # individual draws
                                   bins=bins,          # repeats pre draw
                                   edge=7,             # cutoff for binning at +/- edge
                                   seed=123456789,
                                   steps=50,           # resolution (number of rows until draws)
                                   generators=c("Ziggurat", "MT", "LZLLV",
                                                "GSL", #"V1",
                                                "Gretl", "QL"),
                                   showplot=interactive())
