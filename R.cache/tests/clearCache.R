library("R.cache")

## Use an empty temporary file cache
setCacheRootPath(path=file.path(tempdir()))

## Try to clear it
clearCache(recursive=TRUE, prompt=FALSE)
clearCache(recursive=TRUE, prompt=TRUE)

dirs <- c("tests", "clearCache")
saveCache(1, key=list("clearCache"), dirs=dirs)

clearCache(recursive=FALSE, prompt=TRUE)
clearCache(recursive=TRUE, prompt=TRUE)
clearCache(recursive=FALSE, prompt=FALSE)
clearCache(recursive=TRUE, prompt=FALSE)
