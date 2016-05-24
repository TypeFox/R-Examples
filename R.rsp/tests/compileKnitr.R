library("R.rsp")
library("R.utils") # Arguments

path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.knitr.Rnw")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  if (isCapableOf(R.rsp, "knitr")) {
    outPath <- file.path("LoremIpsum", "knitr.Rnw");
    pathnameR <- compileKnitr(pathname, outPath=outPath, verbose=-10)
    print(pathnameR)
    pathnameR <- Arguments$getReadablePathname(pathnameR)
  }
}

