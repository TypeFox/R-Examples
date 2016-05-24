

cache.dir <- normalizePath(file.path(tempdir(), "harvestr-cache"), mustWork=F)
options(harvestr.cache.dir=cache.dir)
reg.finalizer(emptyenv(), function(...){unlink(cache.dir, TRUE)}, onexit=TRUE)

long_function <- function(wait=15){
  # allow ne to read the ether 
  Sys.sleep(wait)
  paste("Your luck lotto numbers are:", 
    paste(sample(56, 5), collapse=" "),
    '|', sample(46, 1), sep=' ')
}

cache.dir <- getOption("harvestr.cache.dir")
if (file.exists(cache.dir)) unlink(cache.dir, recursive=TRUE, force=TRUE)
t1 <- system.time(a <- withseed(gather(1)[[1]], long_function(20), cache=T))
t2 <- system.time(b <- withseed(gather(1)[[1]], long_function(20), cache=T))

rbind(t1, t2)
identical(a, b)

