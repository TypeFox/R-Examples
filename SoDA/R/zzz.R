
.hasRSPerl <- FALSE

.getMarsData <- function() {
    if(file.exists("mars.csv"))
        eval.parent(quote(mars <- read.csv("mars.csv", skip = 5, as.is = TRUE)))
    else if(file.exists("Examples/mars.csv"))
        eval.parent(quote(mars <- read.csv("Examples/mars.csv", skip = 5, as.is = TRUE)))
    else if(file.exists("Examples/mars.Rdata"))
        eval.parent(quote(load("Examples/mars.Rdata")))
    else
      stop("To run this example, you must download or construct the file \"mars.csv\"; see page 176 in the book")
}

.onLoad <- function(libname, pkgname) {
   ##  .hasRSPerl <<- !inherits(tryCatch(loadNamespace("RSPerl"), error = function(e)e), "error")
   ## if(.hasRSPerl) {
   ##     perlFiles <- readLines(system.file("Perl/perlFiles.txt", package = "SoDA"))
   ##     for(file in perlFiles) {
   ##         ex <- tryCatch(RSPerl::.PerlFile(system.file("Perl", file, package = "SoDA")), error = function(e)e)
   ##         if(inherits(ex, "error"))
   ##           warning("error in running Perl file ", file, "; some Perl subroutines may be missing")
   ##     }
   ## }
}
