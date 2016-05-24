.onAttach <- function(lib, pkg)
{
	packageStartupMessage("
		\n Welcome to enviPick version 1.4 
		\n -> For large file.mzXML, have enough memory allocated to R - e.g. use memory.limit(size=xy) on windows.
		\n -> Ensure to provide centroided & baseline-corrected data.
		\n -> Type webpick() to use the browser UI - use a default browser other than Internet Explorer, e.g. Firefox, Google Chrome 
		\n ");
	
}

