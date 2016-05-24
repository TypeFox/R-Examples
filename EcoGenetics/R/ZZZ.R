#' @keywords internal


.onAttach <- function(...) {
	
	vers <- utils::packageDescription("EcoGenetics", fields = "Version")
  
	textstart<- paste("\n\n",
										"\n", "                ---------------------", "\n",
										"\r                  ","||","EcoGenetics","||",
										"\n", "                ---------------------", "\n",
										"\n", "  Spatial Analysis of Phenotypic, Genotypic and Environmental Data",
										"\n",
										"\n", "  Version",  vers, "\t\t",
										" <leandroroser@ege.fcen.uba.ar>", "\n\n",
										"  Overview: help('EcoGenetics')")
	
	packageStartupMessage(textstart)
}
