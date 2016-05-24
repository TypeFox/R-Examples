#' @useDynLib fmrs

.onAttach <- function(lib, pkg){
  fmrsversion <- utils::packageVersion("fmrs")
  fmrsdate <- utils::packageDescription("fmrs")$Date
  fmrsdescrip <- utils::packageDescription("fmrs")$Description
  fmrsBugReports <- utils::packageDescription("fmrs")$BugReports
  packageStartupMessage(
    paste('fmrs package, Version ',fmrsversion,', Released ',fmrsdate,'\n',
          fmrsdescrip, '\nBugReports: ',fmrsBugReports , sep = "")
  )
}

.onUnload <- function(libpath)
  library.dynam.unload("fmrs", libpath)
