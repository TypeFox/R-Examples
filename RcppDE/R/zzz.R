.onAttach <- function (lib, pkg) {
    packageStartupMessage(paste("Now loading:",
                                "  RcppDE: C++ Implementation of Differential Evolution Optimisation",
                                "  Author: Dirk Eddelbuettel",
                                "Based on:",
                                "  DEoptim (version 2.0-7): Differential Evolution algorithm in R",
                                "  Authors: David Ardia, Katharine Mullen, Brian Peterson and Joshua Ulrich",
                                sep="\n"))
}
