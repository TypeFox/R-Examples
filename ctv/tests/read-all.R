#### Read all CTVs in the package:

library("ctv")

cdir <- system.file("ctv", package="ctv")
ctvs <- list.files(cdir, pattern = "\\.ctv$")
ctvs

rr <- sapply(ctvs,
	     function(ctv) {
		 cat(sprintf("%25s  ",ctv))
		 R <- read.ctv(file.path(cdir, ctv))
		 cat("[ok]\n")
                 R
	     },
	     simplify=FALSE)

for(n in names(rr)) {
    cat(n," :\n", rep.int("=", nchar(n)),"\n", sep='')
    print(rr[[n]])
    cat("--------------------------------------------------------\n")
}
