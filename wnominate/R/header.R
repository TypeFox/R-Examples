.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

   packageStartupMessage("\n## W-NOMINATE Ideal Point Package")
   packageStartupMessage("## Copyright 2006 -", this.year)
   packageStartupMessage("## Keith Poole, Jeffrey Lewis, James Lo, and Royce Carroll")
   packageStartupMessage("## Support provided by the U.S. National Science Foundation")
   packageStartupMessage("## NSF Grant SES-0611974\n")
   #require("pscl", quietly=TRUE)

}

.onUnload <- function(libpath) {
    library.dynam.unload("wnominate", libpath)
}
