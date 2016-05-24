.onAttach <- function (lib, pkg) {
    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION", package="geoR"),
                              fields=c("Title","Version","Date")))

   packageStartupMessage(paste("--------------------------------------------------------------\n",
			       pkg.info["Title"]),"\n",
			 " For an Introduction to geoR go to http://www.leg.ufpr.br/geoR\n",
    paste(" geoR version ", pkg.info["Version"],
              " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""),
			 "--------------------------------------------------------------\n"
    )
}

#".First.lib" <- function(lib, pkg)
#{
#  library.dynam("geoR", package = pkg, lib.loc = lib)
#  messages <- as.logical(ifelse(is.null(getOption("geoR.messages")),
#                                TRUE, getOption("geoR.messages")))
#  if(messages){
#    cat("\n")
#    cat("-------------------------------------------------------------\n")
#    pkg.info <- drop(read.dcf(file=system.file("DESCRIPTION",
#                                package="geoR"),
#                              fields=c("Title","Version","Date")))
#    cat(pkg.info["Title"])
#    cat("\n")
#    cat("For an Introduction to geoR go to http://www.leg.ufpr.br/geoR\n")
#    cat(paste("geoR version ", pkg.info["Version"],
#              " (built on ", pkg.info["Date"], ") is now loaded\n", sep=""))
#    cat("-------------------------------------------------------------\n")
#    cat("\n")
#  }
#  return(invisible(0))
#}

