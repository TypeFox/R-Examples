#' @export
copyRequiredJobFiles.ExperimentRegistry = function(reg1, reg2, id) {
  job = getJob(reg1, id, check.id = FALSE)

  src = getProblemFilePaths(reg1$file.dir, job$prob.id)
  dest = getProblemFilePaths(reg2$file.dir, job$prob.id)
  info("Copying problem files: %s", collapse(src, sep = ", "))
  file.copy(src, dest)

  src = getAlgorithmFilePath(reg1$file.dir, job$algo.id)
  dest = getAlgorithmFilePath(reg2$file.dir, job$algo.id)
  info("Copying algorithm file: %s", src)
  file.copy(src, dest)
}
