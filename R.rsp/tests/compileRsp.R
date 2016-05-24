library("R.rsp")
library("R.utils") # Arguments

# A knitr Rnw file
path <- system.file("rsp_LoremIpsum", package="R.rsp")
pathname <- file.path(path, "LoremIpsum.txt.rsp")
print(pathname)

outPath <- file.path("LoremIpsum", "txt.rsp");
pathnameR <- compileRsp(pathname, outPath=outPath, verbose=-10)
print(pathnameR)
pathnameR <- Arguments$getReadablePathname(pathnameR)

outPath <- file.path("LoremIpsum", "txt.rsp.rsp");
pathnameR2 <- compileRsp(pathname, outPath=outPath, verbose=-10)
print(pathnameR2)
pathnameR2 <- Arguments$getReadablePathname(pathnameR2)
stopifnot(identical(readLines(pathnameR2), readLines(pathnameR)))

outPath <- file.path("LoremIpsum", "txt.rsp.rsp,repeat");
pathnameR3 <- rfile(pathname, workdir=outPath, verbose=-10)
print(pathnameR3)
pathnameR3 <- Arguments$getReadablePathname(pathnameR3)
stopifnot(identical(readLines(pathnameR3), readLines(pathnameR)))
