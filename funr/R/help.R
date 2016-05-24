




generic_help <- function(help_text, script_name = "funr"){
	nm = script_name

	if(missing(help_text))
		help_text = sprintf("
This aims to provide an easy command-line interface to all R functions.

Usage: %s [-h -v] <name of a function> <arguments to the function [<options>=<value>]>

%s -h            Show this help
%s -h <function> Show help for a specific function
%s -v            Show extra verbose prints, for debugging this package
%s <func>        Find and run <function>
%s <func> [args] Run <func> with supplied <arguments>
%s <pkg::func>   Load the package (pkg), and then run <func>

Examples:
    ## Show help for rnorm (random numbers from normal distribution)
    %s -h rnorm
    ## generate 100 random numbers
    %s rnorm n=100

    ## load knitr, then call knit2html to stitch a file
    %s knitr::knit2html <all other arguments>

    ## get an example file from the knitr package
    rmd=$(%s system.file package=knitr fl=examples/knitr-minimal.Rmd)
    ## run knitr on that file
    %s knitr::knit input=$rmd

", nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm, nm)

# 	if(nm == "flowr")
# 		default_help = c(default_help, flow_help())

	return(help_text)
}




#cat("\nFunctions where the arguments are simple objects like numeric/character/logical can be called.",
#		"\nLists become a little more complicated. say x=a,b is converted to a list with elements a,b\n",
#		"\nSome examples:")
#}else{
#cat("\n## Fetch some numbers:\nrfun rnorm n=100",
#			"\n## Fetch files from pacakges:",
#		"\nrmd=$(rfun system.file package=knitr ...=examples/knitr-minimal.Rmd)",
#		"\necho $rmd",
#		"\n## knit this awesome example !",
#		"\nrfun knitr::knit input=$rmd\n")
#}


