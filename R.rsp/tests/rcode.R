library("R.rsp")

# RSP document
pathname <- system.file("exData", "slowcounting.txt.rsp", package="R.rsp")

# Compile to RSP source code script file
metadata <- list(foo="ABC", bar="123")
pathnameR <- rcode(file=pathname, metadata=metadata)
print(pathnameR)
# Assert correct filename
stopifnot(basename(pathnameR) == "slowcounting.txt.R")


# Compile to RSP source code object
metadata <- list(foo="ABC", bar="123")
code <- rcode(file=pathname, metadata=metadata, output=RspSourceCode())
print(code)
