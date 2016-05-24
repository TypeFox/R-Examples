dso.path <- function(fx) {
  # find the path for the dynamic shared objects associated with 
  # the returned object from cxxfunction 
  # 
  # Args:
  #   fx: returned object from cxxfunction in package inline 
  dllinfo <- getDynLib(fx)
  dllinfo[['path']] 
} 

read.dso <- function(path) {
  n <- file.info(path)$size
  readBin(path, what = 'raw', n = n)
} 

cxxfunctionplus <- function(sig = character(), body = character(),
                            plugin = "default", includes = "",
                            settings = getPlugin(plugin), 
                            save.dso = FALSE, ..., verbose = FALSE) {
  fx <- cxxfunction(sig = sig, body = body, plugin = plugin, includes = includes, 
                    settings = settings, ..., verbose = verbose)
  dso.last.path <- dso.path(fx)
  dso.bin <- if (save.dso) read.dso(dso.last.path) else raw(0)
  dso.filename <- sub("\\.[^.]*$", "", basename(dso.last.path)) 
  if (!is.list(sig))  { 
    sig <- list(sig) 
    names(sig) <- dso.filename 
  } 
  dso <- new('cxxdso', sig = sig, dso.saved = save.dso, 
             dso.filename = dso.filename, 
             dso.bin = dso.bin, 
             system = R.version$system, 
             .MISC = new.env()) 
  assign("cxxfun", fx, envir = dso@.MISC)
  assign("dso.last.path", dso.last.path, envir = dso@.MISC)
  return(dso)
} 

# write.dso 
# writeBin(dso, '/tmp/Rtmpdb9w5A/aa.so')
