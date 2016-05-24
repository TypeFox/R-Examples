library("R.rsp")
library("R.utils") # Arguments

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_LoremIpsum")
pathname <- file.path(path, "LoremIpsum.asciidoc.Rnw")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  postprocess <- isCapableOf(R.rsp, "asciidoc");
  outPath <- file.path("LoremIpsum", "asciidoc.Rnw");
  pathnameR <- compileAsciiDocNoweb(pathname, outPath=outPath, postprocess=postprocess, verbose=-10)
  print(pathnameR)
  pathnameR <- Arguments$getReadablePathname(pathnameR)
}

