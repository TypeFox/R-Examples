library("R.rsp")
library("R.utils") # Arguments

# A knitr Rnw file
path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.knitr.Rnw")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  outPath <- file.path("LoremIpsum", "knitr.Rnw-auto");
  pathnameR <- compileRnw(pathname, outPath=outPath, verbose=-10)
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}


# A Sweave Rnw file
path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.Rnw")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  outPath <- file.path("LoremIpsum", "Rnw-auto");
  pathnameR <- compileRnw(pathname, outPath=outPath, verbose=-10)
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}
