as.winsteps <- function(cmd, ifile, pfile, daterun = NULL,
  comptime = NULL) {

  out <- list(cmd = as.wcmd(cmd), ifile = as.ifile(ifile),
    pfile = as.pfile(pfile), daterun = daterun,
    comptime = comptime)

  class(out) <- "winsteps"

  return(out)
}
