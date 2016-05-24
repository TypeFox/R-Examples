library("R.rsp")
library("R.utils") # copyDirectory() and Arguments

path <- system.file(package="R.rsp")
path <- file.path(path, "rsp_LoremIpsum")
pathname <- file.path(path, "LoremIpsum.asciidoc.txt")
print(pathname)

if (Sys.getenv("_R_CHECK_FULL_") != "") {
  if (isCapableOf(R.rsp, "asciidoc")) {
    outPath <- file.path("LoremIpsum", "asciidoc.txt");
    copyDirectory(file.path(path, "figures"), file.path(outPath, "figures"))
    pathnameR <- compileAsciiDoc(pathname, outPath=outPath, verbose=-10)
    print(pathnameR)
    pathnameR <- Arguments$getReadablePathname(pathnameR)
  }
}

