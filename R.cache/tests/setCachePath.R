library("R.cache")
setCachePath <- R.cache:::setCachePath

## Use an empty temporary file cache
tmpdir <- tempdir()
setCacheRootPath(path=tmpdir)

dirs <- c("tests", "readCacheHeader")
setCachePath(dirs=dirs, path=tmpdir)

## Cleanup
clearCache(recursive=TRUE)
