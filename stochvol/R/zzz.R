.onAttach <- function(lib, pkg) {
	if(interactive() || getOption("verbose")){
		packageStartupMessage(sprintf("Package %s (%s) attached. To cite, see citation(\"%s\").", pkg, utils::packageDescription(pkg)$Version, pkg))
	}
}
