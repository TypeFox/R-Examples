#' @method updateRegistry ExperimentRegistry
#' @export
updateRegistry.ExperimentRegistry = function(reg) {
  # update the BatchJobs part first
  updated = NextMethod()
  is.updated = !isFALSE(updated)
  if (is.updated) reg = updated

  version.reg = reg$packages$BatchExperiments$version
  version.pkg = packageVersion("BatchExperiments")
  if (version.reg == version.pkg) {
    return(if (is.updated) reg else FALSE)
  }

  if (version.reg > version.pkg) {
    warningf("The registry has been used with BatchExperiments version %s, installed is version %s. You should update BatchExperiments on this machine.",
             version.reg, version.pkg)
    return(if (is.updated) reg else FALSE)
  }

  # do updates
  if (version.reg < package_version("1.0.767")) {
    path = file.path(reg$file.dir, "problems")
    src = list.files(path, full.names = TRUE, pattern = "_static\\.RData$")
    ids = sub("_static\\.RData$", "", basename(src))
    dest = vapply(ids, function(id) getProblemFilePaths(reg$file.dir, id)["static"], character(1L))
    file.rename(src, dest)

    path = file.path(reg$file.dir, "problems")
    src = list.files(path, full.names = TRUE, pattern = "_dynamic\\.RData$")
    ids = sub("_dynamic\\.RData$", "", basename(src))
    dest = vapply(ids, function(id) getProblemFilePaths(reg$file.dir, id)["dynamic"], character(1L))
    file.rename(src, dest)

    path = file.path(reg$file.dir, "algorithms")
    src = list.files(path, full.names = TRUE)
    ids = sub("\\.RData$", "", basename(src))
    dest = getAlgorithmFilePath(reg$file.dir, ids)
    file.rename(src, dest)
  }

  reg$packages$BatchExperiments$version = version.pkg
  reg
}
