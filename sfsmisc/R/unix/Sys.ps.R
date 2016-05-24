#### Martin Maechler, Aug.2000, originally in /u/maechler/R/MISC/ps.R

### --> ../../man/unix/Sys.ps.Rd to see more comments
### --> ../../man/unix/Sys.ps.Rd to see more comments

## I would really like builtin   Sys.ps() for these

### Sys.ps.cmd() ---> now moved to ../R/

## These only apply to "System V" compatible `ps', not to BSD ones
.Sys.ps.fields <-
  list(POSIX = c("args", "comm", "time", "etime", "nice", "pcpu",
         "pid", "pgid", "ppid", "group", "rgroup", "user", "ruser",
         "tty", "vsz"),

       ## Now the extras, not in above POSIX:
       SunOS = c( "addr", "pri", "c", "rgid", "class", "rss", "f",
         "ruid", "fname", "s", "gid", "sid", "opri", "stime", "osz",
         "uid", "pmem", "wchan"),
       Linux =## These are Linux (RH 6.2):"Docu" at end ..
       c("%cpu", "%mem", "alarm", "blocked", "bsdstart", "bsdtime",
         "c", "caught", "cmd", "command", "cputime", "drs", "dsiz", "egid",
         "egroup", "eip", "esp", "euid", "euser", "f", "fgid", "fgroup",
         "flag", "flags", "fname", "fsgid", "fsgroup", "fsuid", "fsuser",
         "fuid", "fuser", "gid", "ignored", "intpri", "lim", "longtname",
         "lstart", "m_drs", "m_trs", "maj_flt", "majflt", "min_flt", "minflt",
         "ni", "nwchan", "opri", "pagein", "pending", "pgrp", "pmem",
         "pri", "rgid", "rss", "rssize", "rsz", "ruid", "s", "sess", "session",
         "sgi_p", "sgi_rss", "sgid", "sgroup", "sid", "sig", "sig_block",
         "sig_catch", "sig_ignore", "sig_pend", "sigcatch", "sigignore",
         "sigmask", "stackp", "start", "start_stack", "start_time", "stat",
         "state", "stime", "suid", "suser", "svgid", "svgroup", "svuid",
         "svuser", "sz", "timeout", "tmout", "tname", "tpgid", "trs",
         "trss", "tsiz", "tt", "tty4", "tty8", "ucomm", "uid", "uid_hack",
         "uname", "vsize", "wchan")
       )

## Note that  proc.time() gives part of that info better

## command == cmd == args   gives "command + arguments : too long
.Sys.ps.multifields <- c("command", "cmd","args", "lstart")

Sys.ps <-
    function(process = Sys.getpid(),
             fields = c("pid", "pcpu", "time", "vsz", "comm"),
             usefile = length(fields) > 10,
             ps.cmd  = Sys.ps.cmd(),
             verbose = getOption("verbose"),
             warn.multi = verbose || any(fields != "ALL"))
{
    if(!is.null(tp <- attr(ps.cmd,"type")) && tp == "BSD")
        stop("this function cannot work with BSD kind of `ps'.")

    ps.opt <- {
        if(is.numeric(process) && process == round(process))
            paste("-p",process) # PID
        else if(process == "ALL") {
            warning("`process = \"ALL\"' not yet working properly")
            "-e" # all process
          }
        else if(is.character(process) && length(process) == 1)
            paste("-C",process) # Command name
        else stop(paste("invalid `process':",format(process)))
    }
    if(length(ps.opt) > 1)
      warning("Multiple processes : not yet working ...")

    Sys.ps.fields <- c(.Sys.ps.fields $ POSIX,
                       if(any(ii <- Sys.info()["sysname"] ==
                              names(.Sys.ps.fields)))
                       .Sys.ps.fields[ii][[1]])

    if(identical(fields, "ALL"))
      i.field <- TRUE
    else {
      i.field <- pmatch(fields, Sys.ps.fields) # allow abbreviated ones
      if(any(ina <- is.na(i.field))) {
        warning(paste("Dropping invalid field names",
                      fields[ina]))
        i.field <- i.field[!ina]
      }
    }
    fields <- Sys.ps.fields[i.field]
    imult <- !is.na(match(fields, .Sys.ps.multifields))
    if(any(imult) && length(fields) > 1) {
      if(warn.multi)
        warning(paste("Not using `multi fields' ",
                      paste(fields[imult],collapse=",")))
        fields <- fields[!imult]
        imult <- FALSE
    }
    ## Don't use "-w" with  cmd/args, or command : gives space in between
    ## Must use "--width" (GNU ps only) when there are many fields ...
    ## Need temporary file & scan since system cannot get very long
    ## lines ...
    if(usefile)
        ofile <- tempfile("R.Sys.ps")
    cmd <- paste(ps.cmd, ps.opt,
                 "-o", paste(fields, collapse=","),
                 if(usefile) paste(" >", ofile))
    if(verbose) cat("Now calling\n\t",cmd,"\n")
    lines <- system(cmd, intern = !usefile)
    if(usefile) {
        if(lines) warning(paste("system() returned non-0 :",lines))
        lines <- scan(ofile, what = "", sep="\n", quiet = TRUE)## incl header
    }
    if(length(lines) <= 1)
        stop(paste("call returned less than two lines:", lines, sep="\n"))

    r <- sub("^ ","", gsub("[ 	]+"," ", lines))
    ##                     SP & TAB
    if(length(fields) == 1) {
        if(length(r) == 2)
            return(structure(r[2], names = fields))
        else
            warning(paste("Funny result with one `fields': length(r)=",
                          length(r)))
    }
    ## else {
    ll <- strsplit(r, " ")
    d.len <- diff(lenl <- sapply(ll, length))
    if(lenl[1] == length(fields))
        ## use fields!
        ll[[1]] <- fields
    else
        warning(paste("Number returned headers =", lenl[1], " != ",
                      "#{fields} =", length(fields)))
    if(d.len) { # names and result differ
        warning(paste("Lengths differ:",
                      paste(lenl, collapse=",")))
    }
    r <- c(ll[[2]], rep(NA, max(0,-d.len)))
    names(r) <-
        if( d.len > 0) c(ll[[1]], rep(".x.",d.len)) else ll[[1]][1:lenl[2]]
    r
    ##}
}

Sys.sizes <- function(process = Sys.getpid(),
                      ps.cmd  = Sys.ps.cmd())
{
  ## For both Solaris and GNU(Linux);  GNU/Linux additionally has dsize

  if(!is.null(tp <- attr(ps.cmd,"type")) && tp == "BSD") {
    ## a *real* hack [needed for Linux 2.0 or SunOS 4.x ..]
    r <- system(paste(ps.cmd,"m",process), intern = TRUE)[1:2] # 2 lines
    r <- strsplit(r,"  *")
    hd <- r[[1]]; hd <- hd[hd != "" & hd != "COMMAND"]
    i <- match(c("RSS","DRS"), hd)
    r <- structure(r[[2]][i], names = hd[i])
  }
  else { ## proper "System V like"  ps :
    r <- Sys.ps(process, c("rss","vsz"))
  }
  storage.mode(r) <- "integer"
  r
}

if(Sys.info()[["sysname"]] == "Linux") { ##----- Linux-only ----

 Sys.procinfo <- function(procfile)
 {
   l2 <- strsplit(readLines(procfile),"[ \t]*:[ \t]*")
   r <- sapply(l2[sapply(l2, length) == 2],
               function(c2)structure(c2[2], names= c2[1]))
   attr(r,"Name") <- procfile
   names(r) <- make.names(names(r), unique = TRUE) # <- so the result can be name-indexed!
   class(r) <- "simple.list"
   r
 }

} else { ## non-Linux "unix" -- including MacOS X "Darwin"

 Sys.procinfo <- function(procfile) {
     stop("Sys.procinfo() is not yet implemented for non-Linux unix-alikes")
 }

}

Sys.cpuinfo <- function() Sys.procinfo("/proc/cpuinfo")
Sys.meminfo <- function() Sys.procinfo("/proc/meminfo")

Sys.MIPS <- function() as.numeric(Sys.cpuinfo()["bogomips"])

Sys.memGB <- function(kind = "MemTotal") {
    mm <- drop(read.dcf("/proc/meminfo", fields=kind))
    if(any(is.na(mm))) stop("Non-existing 'kind': ", names(mm)[is.na(mm)][1])
    if(!all(grepl(" kB$", mm)))
        stop("Memory info ", dQuote(kind), " is not returned in 'kB' aka kiloBytes")
    ## return memory in giga bytes
    as.numeric(sub(" kB$", "", mm)) / (1000 * 1024)
}
