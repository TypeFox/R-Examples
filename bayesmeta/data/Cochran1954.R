#
#  W.G. Cochran. The combination of estimates from different experiments.
#  Biometrics, 10(1):101-129, 1954.
#
#  (Example 3, p.119)
#

Cochran1954 <- data.frame("observer"   = as.character(1:7),
                          cbind("mean" = c(183.2, 149.0, 154.0, 167.2, 187.2, 158.0, 143.0),
                                "se2"  = c(117.0, 8.1, 235.9, 295.0, 1064.6, 51.2, 134.0)),
                          stringsAsFactors=FALSE)
