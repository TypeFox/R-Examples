winsteps <- function(cmd, cmdfile = "cmdfile", outfile = "outfile",
  ifile = "ifile", pfile = "pfile", newdir = getwd(),
  run = TRUE, windir = "winsteps") {

  olddir <- getwd()
  setwd(newdir)

  if(run) {
    if(!missing(cmd))
      write.wcmd(cmd, filename = cmdfile)

    systemcommand <- paste(windir, "BATCH=YES", cmdfile, outfile,
      paste("PFILE=", pfile, sep = ""),
      paste("IFILE=", ifile, sep = ""))
    gc(FALSE)
    time1 <- proc.time()
    outval <- system(systemcommand)
    time2 <- proc.time()
    if(outval != 0)
      stop("Winsteps not run - error sending command file")
    else
      cat("\nCommand file sent to Winsteps\n\n")
  }

  out <- as.winsteps(cmd = read.wcmd(cmdfile),
    ifile = read.ifile(ifile), pfile = read.pfile(pfile),
    daterun = date(), comptime = time2 - time1)

  if(cmdfile == "cmdfile")
    unlink("cmd")
  if(pfile == "pfile")
    unlink("pfile")
  if(ifile == "ifile")
    unlink("ifile")
  if(outfile == "outfile")
    unlink(outfile)

  setwd(olddir)

  return(out)
}
