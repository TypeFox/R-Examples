library("R.cache")

## Use an empty temporary file cache
setCacheRootPath(path=file.path(tempdir()))

dirs <- c("tests", "readCacheHeader")

for (compress in c(FALSE, TRUE)) {
  pathname <- saveCache(1, key=list("readCacheHeader"), dirs=dirs, compress=compress)

  for (byName in c(FALSE, TRUE)) {
    if (byName) {
      hdr <- readCacheHeader(pathname)
    } else {
      con <- gzfile(pathname, open="rb")
      hdr <- readCacheHeader(con)
      close(con)
    }
    str(list(pathname=pathname, hdr=hdr))
  }
}

## Cleanup
clearCache(recursive=TRUE)
