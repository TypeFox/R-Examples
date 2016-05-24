.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

   packageStartupMessage("\n## alpha-NOMINATE Ideal Point Package")
   packageStartupMessage("## Copyright 2013 - ", this.year)
   packageStartupMessage("## Royce Carroll, Christopher Hare, Jeffrey B. Lewis, James Lo, Keith Poole, and Howard Rosenthal")
   packageStartupMessage("## Support provided by the U.S. National Science Foundation")
   packageStartupMessage("## NSF Grant SES-0611974\n")

}

.onUnload <- function(libpath) {
    library.dynam.unload("anominate", libpath)
}
