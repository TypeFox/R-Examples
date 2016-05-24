.onAttach <- function(...) {
  mylib <- dirname(system.file(package = "anchors"))
  ver <- packageDescription("anchors", lib.loc = mylib)$Version
  builddate <- packageDescription("anchors", lib.loc = mylib)$Date
  packageStartupMessage("\n##  anchors (Version ", ver,", Build Date: ", builddate, ")\n",
	"##  See http://wand.stanford.edu/anchors for additional documentation and support.\n\n")
#  cat("## Contact Jonathan Wand <wand(at)stanford.edu> with comments about anchors.\n", sep="")
#  cat("##\n", sep="")
}
