#
#  G.W. Snedecor and W.G. Cochran. Statistical Methods. Iowa State University Press,
#  Ames, IA, USA, 6th edition, 1967.
#
#  J. Hartung, G. Knapp, and B.K. Sinha. Statistical meta-analysis with applications.
#  Wiley, Hoboken, NJ, USA, 2008.
#
#  (example 10.18.1, p.290, in Snedecor/Cochran,
#   table 7.3, p.90, in Hartung/Knapp/Sinha)
#

SnedecorCochran <- data.frame("no"   = paste("no.", 1:6, sep=""),
                              cbind("n"    = c(5,2,7,5,7,9),
                                    "mean" = c(41.2,  64.5,  56.3,   39.6,  67.1,  53.2),
                                    "var"  = c(35.14, 30.25, 18.99, 101.06, 38.64, 27.72)),
                              stringsAsFactors=FALSE)
