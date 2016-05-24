.noGenerics <- TRUE

.onAttach <- function(lib, pkg){
  packageStartupMessage(paste("Now loading:",
                                "  bbo: R Implementation of Biogeography-Based Optimization",
                                "  Author: Sarvesh Nikumbh",
                                "Based on:",
                                "  D. Simon, \u0022Biogeography-Based Optimization,\u0022 IEEE Transactions on Evolutionary Computation, vol. 12, no. 6, pp. 702-713, December 2008",
                                sep="\n")
			 )
}
