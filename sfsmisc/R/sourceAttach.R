##' @title "Source + Attach" an R source file
##' @author Martin Maechler, 29 Jul 2011
sourceAttach <- function(file, pos = 2,
			 name = paste(abbreviate(gsub(fsep, "", dirname(file)), 12,
						 method="both.sides"),
				      basename(file), sep=fsep),
			 keep.source = getOption("keep.source.pkgs"),
			 warn.conflicts = TRUE)
{
    ENV <- new.env()
    sys.source(file, envir = ENV, keep.source = keep.source)# also checks file
    fsep <- .Platform$file.sep # for default 'name' :
    attach(ENV, pos=pos, name=name, warn.conflicts=warn.conflicts)
}
