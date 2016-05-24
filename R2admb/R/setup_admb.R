#'Set up AD Model Builder environment variables
#'
#'Attempts to set environment variables so that AD Model Builder will "just
#'work" when run from inside R
#'
#'(1) If the environment variable ADMB_HOME is not already set and admb_home is
#'not specified, this function will try to set it sensibly. (I.e., on Unix
#'systems, it will run a "locate" command (if one is available) to try to find
#'the binaries, and thereafter check if they are installed in the default
#'location (/usr/local/admb); on Windows it will assume they are installed in
#'the default location (C:/ADMB).) (2) If ADMB_HOME is set and admb_home is not
#'specified, it will leave the original setting alone. (3) If admb_home is
#'specified, it will set the environment variable ADMB_HOME to that value.
#'
#'The function also prepends the admb_home value to the PATH variable.
#'
#'@usage  setup_admb(admb_home)
#' 
#'        clean_admb(fn,which=c("sys","output"))
#'@aliases setup_admb clean_admb
#'@param admb_home (character) directory containing AD Model Builder binary
#'files
#'@param fn (character) base name of ADMB model files
#'@param which what to remove: any combination of "sys" (system), "input",
#'"output", or "all" or "none"
#'@return A character vector containing the name of the current ADMB_HOME.
#'@author Ben Bolker
#'@export clean_admb setup_admb
#'@keywords misc environment
#'@examples
#'
#'  orig <- Sys.getenv("ADMB_HOME")
#' ## this doesn't make sense but won't break anything
#' ## until you actually try to run AD Model Builder
#'  setup_admb("elsewhere")   
#'  Sys.setenv(ADMB_HOME="") ## erase environment variable
#'\dontrun{
#'  setup_admb()              ## auto-locate (fails if ADMB not found)
#'}
#'  Sys.setenv(ADMB_HOME=orig) ## restore sanity
#'
setup_admb <- function(admb_home) {
    ## check whether already set up
    sys_home <- Sys.getenv("ADMB_HOME")
    if (missing(admb_home)) { ## not passed to function
        if (nchar(sys_home)>0) {
            ## previously defined
            admb_home <- sys_home
        } else {
            ## not defined & not passed: try to guess
            admb_home <- ""
            if (.Platform$OS.type=="unix") {
                ## unix.  Should (1) check that locate really exists,
                ## (2) check that it finds something, (3) try to
                ## use 'default' location??
                admb_home <- suppressWarnings(system("locate bin/admb | grep bin/admb$",intern=TRUE))
                if (length(admb_home)>0) {
                    admb_home <- gsub("/bin/admb$","",admb_home)
                } else if (file.exists("/usr/local/admb")) {
                    ## try default location
                    admb_home <- "/usr/local/admb"
                } else if (file.exists("/Applications/ADMBTerminal.app/admb/admb")) {
                    ## try default location for MacOSx applicationelse {
                    admb_home <- "/Applications/ADMBTerminal.app/admb"
                } else
                    admb_home <- ""
                ## n.b. extra slash at end of ADMB_HOME is **VERY BAD** **VERY CONFUSING**
                ##  provokes weird behavior where "bin/sedd..." turns into "binsedd..." ???
                if (length(admb_home)>1) {
                    warning("'locate' found more than one instance of bin/admb: using last")
                    ## FIXME: query user for which one to use?
                    admb_home <- admb_home[length(admb_home)]
                }
            }
            if (.Platform$OS.type=="windows") {
                ## default location from IDE setup
                admb_home <- "c:/admb"
            }
            if (nchar(admb_home)==0) stop("couldn't guess ADMB_HOME location,",
                     "you will have to configure it manually")
        }
    }
    Sys.setenv(ADMB_HOME=admb_home)
    path=Sys.getenv("PATH")
    if (.Platform$OS.type=="windows") {
        pathsepchr <-  ";"
        Sys.setenv(PATH=paste(
                   paste(paste(admb_home,c("bin","utilities"),sep="/"),
                         collapse=pathsepchr),
                   "C:/MinGW/bin",
                   path,
                   sep=pathsepchr))
        
        
    } else {
        pathsepchr <-  ":"
        Sys.setenv(PATH=paste(paste(admb_home,"bin",sep="/"),
                   path,
                   sep=pathsepchr))
    }
    return(admb_home)
}

clean_admb <- function(fn,which=c("sys","output")) {
    if (length(which)==1) {
        if (which=="none") return()
        if (which=="all") which <- c("sys","input","output","gen")
    }
    
    sys.ext <- c("bar","bgs","cpp","ecm","eva","htp","luu","mc2","mcm","o","rep","rhes",
                 "luu","mc2","mcm","tpl.bak","out","cout","shess")
    sys.files <- paste(fn,sys.ext,sep=".")
    gen.files <- list.files(pattern="_gen(\\.tpl)*")
    sys.other <- c("eigv.rpt","fmin.log","variance","sims",
                   "hesscheck","hessian.bin","dgs2","diags",
                   paste("admodel",c("dep","hes","cov"),sep="."),
                   list.files(pattern="xx.*.tmp"),
                   list.files(pattern=".*f1b2list.*"),
                   list.files(pattern=paste(fn,"\\.[bpr][0-9]+",sep="")))
    ## FIXME: clean up abandoned 'buffer' files too
    ## f1b2list etc.
    input.ext <- c("pin","dat")
    input.files <- paste(fn,input.ext,sep=".")
    output.ext <- c("log","cor","std", "par","psv","hst","prf","mcinfo")
    output.files <- paste(fn,output.ext,sep=".")
    output.files <- c(output.files,list.files(pattern="\\.plt$"))
    if ("sys" %in% which) unlink(c(sys.files,sys.other))
    if ("input" %in% which) unlink(input.files)
    if ("output" %in% which) unlink(output.files)
    if ("gen" %in% which) unlink(gen.files)
}

## @rm -vf $(PROGRAM_NAME){.htp,.cpp,.std,.rep,.b??,.p??,.r??,.cor,.eva,.log,.rhes,.luu,}
## @rm -vf admodel{.cov,.dep,.hes} eigv.rpt fmin.log variance hess*

