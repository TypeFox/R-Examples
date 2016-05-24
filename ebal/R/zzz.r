.onAttach <-
function (libname, pkgname)
{
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
   packageStartupMessage("##\n## ebal Package: Implements Entropy Balancing.\n")
   packageStartupMessage("## See http://www.stanford.edu/~jhain/ for additional information.\n\n")

}

