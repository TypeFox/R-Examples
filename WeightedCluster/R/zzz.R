
.onAttach <- function(libname, pkgname){
	suppressWarnings(descr <- utils::packageDescription("WeightedCluster"))
	if(utils::packageVersion("WeightedCluster")$minor %% 2 == 0) {
		state <- "stable"
	}
	else {
		state <- "development"
	}
	if(!is.null(descr$Built)){
		builtDate <- paste(" (Built: ", strsplit(strsplit(descr$Built, ";")[[1]][3], " ")[[1]][2], ")", sep="")
	}else{
		builtDate <- ""
	}
	packageStartupMessage("This is WeightedCluster ", state, " version ", descr$Version, builtDate)
	packageStartupMessage('\nTo get the manuals, please run:')
	packageStartupMessage('   vignette("WeightedCluster") ## Complete manual in English')
	packageStartupMessage('   vignette("WeightedClusterFR") ## Complete manual in French')
	packageStartupMessage('   vignette("WeightedClusterPreview") ## Short preview in English')
	packageStartupMessage("\nTo cite WeightedCluster in publications please use:")
	packageStartupMessage("Studer, Matthias (2013). WeightedCluster Library Manual: A practical guide to")
	packageStartupMessage("   creating typologies of trajectories in the social sciences with R.")
	packageStartupMessage("   LIVES Working Papers, 24. doi: 10.12682/lives.2296-1658.2013.24")
	
}
