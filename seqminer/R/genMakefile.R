#' Create a new job
#'
#' @param id character, job ids.
#' @param cmd character, commands to run
#' @param outFile character, the output file names after command are run successfully
#' @param depend character vector, specify the prerequisite files (e.g. outFile from other jobs)
#' @export
#' @examples
#' j1 <- newJob('id1', 'cmd out1', 'out1')
#' j2 <- newJob('id2', 'cmd out2', 'out2', depend = 'out1')
newJob <- function(id, cmd, outFile, depend = NULL) {
  if (!all(grepl(pattern = outFile, x= cmd))) {
    warning(gettextf("Cannot locate outFile in cmd"))
  }
  outFile = paste(outFile, collapse = " ")
  depend  = paste(depend, collapse = " " )
  ret <- list(id = id, cmd = cmd, outFile = outFile, dep = depend)
  class(ret) <- "job"
  ret
}

#' Create a new workflow
#'
#' @param name character, specify the name of the workflow
#' @export
#' @examples
#' w <- newWorkflow("wf")
newWorkflow <- function(name) {
  ret <- list(job = list(), name = name, job.id = character())
  class(ret) <- "workflow"
  ret
}

#' Add a job to a workflow
#'
#' @param wf a variable of workflow class
#' @param job a variable of job class
#' @export
#' @examples
#' j1 <- newJob('id1', 'cmd out1', 'out1')
#' j2 <- newJob('id2', 'cmd out2', 'out2', depend = 'out1')
#' w <- newWorkflow("wf")
#' w <- addJob(w, j1)
#' w <- addJob(w, j2)
#' writeWorkflow(w, "Makefile")
addJob <- function(wf, job) {
  if (!inherits(job, "job")) {
    stop("Only job class can be added")
  }
  if (job$id %in% wf$job.id ) {
    stop("Duplicated job id detected")
  }
  wf$job[[length(wf$job)+1]] <- job
  wf$job.id <- c(wf$job.id, job$id)
  wf
}

#' Export workflow to Makefile
#'
#' @param wf a variable workflow class
#' @param outFile character, typically named "Makefile"
#' @export
#' @examples
#' j1 <- newJob('id1', 'cmd out1', 'out1')
#' j2 <- newJob('id2', 'cmd out2', 'out2', depend = 'out1')
#' w <- newWorkflow("wf")
#' w <- addJob(w, j1)
#' w <- addJob(w, j2)
#' writeWorkflow(w, "Makefile")
writeWorkflow <- function(wf, outFile) {
  if (!inherits(wf, "workflow")) {
    stop("Wrong input")
  }
  ## make sure deps are a subset of outfile
  deps <- do.call(c, lapply(wf$job, function(x) {strsplit(x$dep, " ")[[1]]}))
  outFiles <- do.call(c, lapply(wf$job, function(x) {strsplit(x$outFile, " ")[[1]]}))
  allJobs <- do.call(c, lapply(wf$job, function(x) {x$id}))
  nonMetDeps <- setdiff(deps, outFiles)
  if (length(nonMetDeps) > 0) {
    message(gettextf("Workflow depends on these files: %s", paste(nonMetDeps, collapse = ", ")))
  }

  cat(".DELETE_ON_ERROR:\n", file = outFile, append = FALSE)
  cat("all: ", paste(outFiles, collapse = " "), "\n", file = outFile, append = TRUE)
  for (l in wf$job) {
    # print(l)
    cat("## ", l$id, "\n", file = outFile, append = TRUE)
    cat(l$outFile, " : ", l$dep, "\n",  file = outFile, append = TRUE )
    for (cmd1 in l$cmd) {
      cat("\t", l$cmd, "\n", file = outFile, append = TRUE)
    }
  }
  message(gettextf("Makefile generated: %s", outFile))
  message("Run using: make -f ", normalizePath(outFile))
}

if (FALSE) {
  w <- newWorkflow("wf")
  inFile <- "input"
  nSplit <- 10
  split.out.file <- sprintf('out%d.finish', seq(nSplit))
  for (i in 1:nSplit) {
    j <- newJob(id = sprintf('id%d', i),
                cmd = sprintf('Rscript ... out%d.finish', i),
                outFile = split.out.file[i],
                depend = inFile)
    w <- addJob(w, j)
  }
  j <- newJob(id = 'combine',
              cmd = 'Rscript ... out1 out2 ... out2000',
              outFile = 'output',
              depend = split.out.file)
  w <- addJob(w, j)
  writeWorkflow(w, "tmp.mf")
}
