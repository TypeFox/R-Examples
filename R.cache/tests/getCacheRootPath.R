library("R.cache")
setupCacheRootPath <- R.cache:::setupCacheRootPath

oopts <- options()

tmpdir <- tempdir()
setCacheRootPath(path=tmpdir)

print(getCacheRootPath())

setupCacheRootPath(defaultPath=tmpdir)
path <- getCacheRootPath(NULL)
print(path)
stopifnot(identical(path, tmpdir))

options("R.cache::rootPath"=NULL)
print(getCacheRootPath())

options("R.cache.path"="foo")
print(getCacheRootPath())

options("R.cache.path"=NULL)
print(getCacheRootPath())

oenv <- Sys.getenv("R_CACHE_PATH")
Sys.setenv("R_CACHE_PATH"="")
print(getCacheRootPath())

path <- getCacheRootPath(NULL)
print(path)
stopifnot(is.null(path))
setupCacheRootPath(defaultPath=tmpdir)

path <- getCacheRootPath(NULL)
print(path)
stopifnot(identical(path, tmpdir))


## Cleanup
options(oopts)
args <- list(oenv)
names(args) <- "R_CACHE_PATH"
do.call(Sys.setenv, args)
