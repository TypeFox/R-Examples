library("R.rsp")
library("R.utils") # Arguments

path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.Rnw")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  outPath <- file.path("LoremIpsum", "Rnw");
  pathnameR <- compileSweave(pathname, outPath=outPath, verbose=-10)
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}
