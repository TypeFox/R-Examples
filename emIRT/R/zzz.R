.onAttach <- function(...) {

   mydate <- date()
   x <- regexpr("[0-9]{4}", mydate)
   this.year <- substr(mydate, x[1], x[1] + attr(x, "match.length") - 1)

  packageStartupMessage("\n## emIRT: Fast Estimation of Ideal Points with Massive Data")
   packageStartupMessage("## 2014 - ", this.year)
   packageStartupMessage("## Kosuke Imai, James Lo, and Jonathan Olmsted\n\n")
}

.onUnload <- function(libpath) {
    library.dynam.unload("emIRT", libpath)
}
