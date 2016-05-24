.onAttach <- function(libname, pkgname) {
  pkg.info <- utils::packageDescription('sdm') 
  packageStartupMessage(paste("sdm ", pkg.info[["Version"]], " (", pkg.info["Date"], ")", sep=""))
  #.containers_env <<- new.env()
  #assign('a',.replicateMethods$new(),envir = .containers_env)
  if (.is.installed('dismo')) {
    jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
    .sdmOptions$addOption('maxJar',file.exists(jar))
  }
  .sdmOptions$addOption('sdmLoaded',FALSE)
  invisible(0)
}

.onUnload <- function(libpath) {
  if (".sdmMethods$userFunctions" %in% search()) detach('.sdmMethods$userFunctions')
  invisible(0)
}

