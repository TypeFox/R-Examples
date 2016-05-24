.restartCalibration =  function(control, type="restart") {
  res.file = paste0(control$restart.file, ".", type)
  return(file.exists(res.file))
}

.getRestart = function(control, ...) {
  res.file = paste0(control$restart.file, ".restart")
  load(res.file)
  if(!exists("opt")) stop("Restart file ", res.file, " is not appropiate.")
  if(!exists("trace")) stop("Restart file ", res.file, " is not appropiate.")
  if(!any(class(opt)=="optimES.restart")) stop("Restart file ", res.file, " is not appropiate.")
  opt   = get("opt")
  trace = get("trace")
  return(list(opt=opt, trace=trace))
}

.getResults = function(control, ...) {
  res.file = paste0(control$restart.file, ".results")
  load(res.file)
  if(!exists("output")) stop("Restart file ", res.file, " is not appropiate.")
  if(!any(class(output)=="calibrar.results")) stop("Restart file ", res.file, " is not appropiate.")
  class(output) = "list"
  return(output)
}

.createRestartFile = function(opt, trace, control) {
  if(is.null(control$restart.file)) return(invisible())
  if((opt$gen%%control$REPORT)!=0) return(invisible())
  res.file = paste0(control$restart.file, ".restart")
  class(opt) = c("optimES.restart", class(opt))
  force(trace)
  save(list=c("opt", "trace"), file = res.file, compress=FALSE)
  return(invisible())
}

.createOutputFile = function(output, control) {
  if(is.null(control$restart.file)) return(invisible())
  res.file = paste0(control$restart.file, ".results")
  class(output) = c("calibrar.results", class(output))
  save(output, file = res.file, compress=FALSE)
  suppressWarnings(file.remove(paste0(control$restart.file, ".restart")))
  return(invisible())
}
