.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0( "Warning: This package was not developed by authors affiliated with the Institute for Health Metrics and Evaluation. This is an open source replication of the Tariff algorithm. Unintentional discrepancies may exist.\n",

		 "To cite the Tariff package, please use this function:\n",
		 "citation(package = \"Tariff\")\n",

      	 "\nTo cite the original Tariff method, please use:\n",
         "James, S. L. and A. D. Flaxman and C. J. Murray and Consortium Population Health Metrics Research (2011). ",
         "Performance of the Tariff Method: validation of a simple additive algorithm for analysis of verbal autopsies ",
         "Population Health Metrics, 9(1), pp.1-16.\n"),
      domain = NULL,  appendLF = TRUE )
}
