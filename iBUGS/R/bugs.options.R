#' Set or Query Options to Run WinBUGS/OpenBUGS/JAGS
#'
#' Set or query options in the bugs/jags function.
#'
#'
#' @param \dots any options can be defined, using \sQuote{name = value} or by
#' passing a list of such tagged values. See the examples below.
#' @return A list containing all the options when the query parameter was
#' empty; or the option value when a single parameter was queried. The old
#' option list will be returned if any new options were set.
#' @note The only option that is not in \code{\link[R2WinBUGS]{bugs}} is
#' \sQuote{model.name} which specifies the R object name to be used to store
#' the returned values from \code{\link[R2WinBUGS]{bugs}}. There are several
#' options in \code{\link[R2WinBUGS]{bugs}} and \code{\link[R2jags]{jags}} are
#' omitted since they are not very usefull.
#' @author Yihui Xie <\url{http://yihui.name}>
#' @seealso \code{\link[R2WinBUGS]{bugs}}, \code{\link[R2jags]{jags}}
#' @keywords misc
#' @export
#' @examples
#'
#' \dontrun{
#'
#' ## list all the options
#' bugs.options()
#'
#' ## set options
#' bugs.options(n.chains = 1, n.iter = 10000)
#'
#' ## all these can be set via the GUI: iBUGS()
#'
#' }
#'
bugs.options = function(...) {
  if (is.null(getOption("iBUGS")) || !length(grep(getOption("iBUGS")$program, getOption("iBUGS")$bugs.directory))) {
    bugs.directory = ""
    program = ifelse(is.null(getOption("iBUGS")), "", getOption("iBUGS")$program)
    if (!is.null(getOption("iBUGS"))) {
      ## looking for (Win|Open)BUGS|JAGS
      if (.Platform$OS.type == "windows") {
        if (nzchar(prog <- Sys.getenv("ProgramFiles")) && length(bugs.dir <- list.files(prog, "^((Open|Win)BUGS|JAGS).*")) && length(bugs.exe <- dirname(list.files(file.path(prog, 
                                                                                                                                                                              bugs.dir), pattern = "((Open|Win)BUGS|jags).*\\.exe$", full.names = TRUE, recursive = TRUE)))) {
          if (length(grep(program, bugs.dir))) {
            bugs.directory = bugs.exe[grep(program, bugs.exe)][length(grep(program, bugs.exe))]
            ## OpenBUGS and WinBUGS can be chosen by users under Windows; ignore multiple directories if BUGS installed in multiple places,
            ## choose the latest one
            if (program == "JAGS") 
              require(R2jags)
          } else gmessage(paste(program, " is not installed. Please choose ", bugs.dir, " for Program!"))
        } else gmessage("OpenBUGS, WinBUGS and JAGS are not found installed.")
      } else {
        if (!require(R2jags)) 
          gmessage("Please install JAGS software and R2jags package.")
        program = bugs.directory = "JAGS"
        # I'm considering making R2jags the 'Depends' of iBUGS.
      }
    }
    data = unlist(sapply(grep("^[^(package:)]", search(), value = TRUE), ls))
    inits = NULL
    parameters.to.save = ""
    model.file = "model.bug"
    n.chains = 3
    n.iter = 2000
    n.burnin = floor(n.iter/2)
    n.sims = 1000
    n.thin = max(1, floor(n.chains * (n.iter - n.burnin)/n.sims))
    bin = as.integer((n.iter - n.burnin)/n.thin)
    debug = FALSE
    DIC = TRUE
    digits = 5
    codaPkg = FALSE
    working.directory = NULL
    clearWD = FALSE
    bugs.seed = NULL
    summary.only = FALSE
    save.history = !summary.only
    over.relax = FALSE
    model.name = "bugs.model"
    # Arguments for JAGS only
    jags.seed = 123
    refresh = n.iter/50
    progress.bar = "text"
    mf = list(data = data, inits = inits, parameters.to.save = parameters.to.save, model.file = model.file, n.chains = n.chains, n.iter = n.iter, 
              n.burnin = n.burnin, n.thin = n.thin, n.sims = n.sims, bin = bin, debug = debug, DIC = DIC, digits = digits, codaPkg = codaPkg, 
              bugs.directory = bugs.directory, program = program, working.directory = working.directory, clearWD = clearWD, bugs.seed = bugs.seed, 
              summary.only = summary.only, save.history = save.history, over.relax = over.relax, jags.seed = jags.seed, refresh = refresh, 
              progress.bar = progress.bar, model.name = model.name)
    options(iBUGS = mf)
  } else mf = getOption("iBUGS")
  lst = list(...)
  if (length(lst)) {
    if (is.null(names(lst)) & !is.list(lst[[1]])) {
      getOption("iBUGS")[unlist(lst)][[1]]
    } else {
      omf = mf
      mc = list(...)
      if (is.list(mc[[1]])) 
        mc = mc[[1]]
      if (length(mc) > 0) 
        mf[pmatch(names(mc), names(mf))] = mc
      options(iBUGS = mf)
      invisible(omf)
    }
  } else {
    getOption("iBUGS")
  }
}
