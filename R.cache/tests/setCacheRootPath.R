library("R.cache")

tmpdir <- tempdir()
setCacheRootPath(path=tmpdir)
setCacheRootPath(path=file.path(tmpdir, "subdir"))
setCacheRootPath(path=tmpdir)

## Cleanup
clearCache(recursive=TRUE)
