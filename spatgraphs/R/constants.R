#' Supported graphs constants
#'
#'
SG_GRAPH_PARAMETERS <- list(geometric = list(R="numeric>0"),
                            knn = list(k="integer>0"),
                            mass_geometric=list(mass="numeric vector of sizes"),
                            markcross=list(mass="numeric vector of sizes"),
                            gabriel=list(),
                            MST=list(),
                            SIG=list(),
                            RST=list(center="coordinates of the center"),
                            RNG=list(),
                            CCC=list(types="factor vector of types")
                            )
