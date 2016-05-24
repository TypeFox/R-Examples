library(gtools)

## Examples from man page
Treatment <- c("Control", "Asprin 10mg/day", "Asprin 50mg/day",
               "Asprin 100mg/day", "Acetomycin 100mg/day",
               "Acetomycin 1000mg/day")

stopifnot( mixedorder(Treatment)==c(5, 6, 2, 3, 4, 1) )


x <- rev(c("AA 0.50 ml", "AA 1.5 ml", "AA 500 ml", "AA 1500 ml",
           "EXP 1", "AA 1e3 ml", "A A A", "1 2 3 A", "NA", NA, "1e2",
           "", "-", "1A", "1 A", "100", "100A", "Inf"))

stopifnot( mixedorder(x)==c(7, 11, 4, 5, 3, 8, 2, 1, 6, 12, 18, 17, 16, 13, 15, 14, 10, 9) )

## Bug reported by Aaron Taudt on 2014-03-01

tmp <- c("uniresult_simulated_H3k27ac_binsize_200_chr1.RData",
         "uniresult_simulated_H3k27me3_binsize_200_chr1.RData",
         "uniresult_simulated_H3k36me3_binsize_200_chr1.RData",
         "uniresult_simulated_H3k4me3_binsize_200_chr1.RData",
         "uniresult_simulated_H3k9me3_binsize_200_chr1.RData")

stopifnot( mixedorder(tmp)==c(4, 5, 1, 2, 3) )


