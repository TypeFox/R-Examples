#' @title Plots Y traces of multiple optimization paths
#' 
#' @description Can be used for only single-objective optimization paths.
#' Useful to compare runs of different algorithms on the same optimization problem.
#' You can add your own ggplot layers to the resulting plot object.
#'
#' @param opt.paths [\code{list}]\cr
#'   List of \code{OptPath} objects
#' @param over.time [\code{character}]\cr
#'   Should the traces be plotted versus the iteration number or the cumulated
#'   execution time? For the later, the opt.path has to contain a extra coloum
#'   names exec.time. Possible values are dob and exec.time, default is \code{dob}.
#' @return 
#'   ggplot2 plot object

renderYTraces = function(opt.paths, over.time = "dob") {
  
  assertList(opt.paths, types = "OptPath", min.len = 1L)
  assertChoice(over.time, choices = c("dob", "exec.time"))
  run.name = names(opt.paths)
  y.name = NULL
  minimize = NULL
  data = data.frame()
  # combine y traces for this algo + add algo / repl + do some sanity checks
  fronts = lapply(seq_along(opt.paths), function(j) {
    run = opt.paths[[j]]
    name = names(as.data.frame(run, include.x = FALSE, include.rest = FALSE))
    if (j == 1L) {
      y.name <<- name
      minimize <<- run$minimize
    }
    if (length(name) != 1L) {
      stopf("Must always have 1 objective in opt path. But found: %i", length(name))
    }
    if (!y.name == name) {
      stopf("Must always have the same objective in opt path: %s (first ones taken). But found here: %s",
        y.name, name)
    }
    if (!minimize == run$minimize) {
      stopf("Must always have the same 'minimize' settings for objective in opt path: %s (first one taken).
          But found here: %s", minimize, run$minimize)
    }
    if (over.time == "dob")
      df = data.frame(
        y = getOptPathY(op = run),
        time = getOptPathDOB(op = run)
      )
    if (over.time == "exec.time") {
      times = seq(0, sum(getOptPathExecTimes(run)), length.out = 128)
      df = getOptPathColAtTimes(run, times)
      df = df[, c(y.name, "time")]
    }

    # add column for algorithm and replication
    cbind(df, .algo = run.name[j])
  })
  
  data = do.call(rbind, fronts)
  mean.data = aggregate(data[[y.name]], by = list(data$time, data$.algo), FUN = mean)
  names(mean.data) = c("time", ".algo", y.name)
  
  pl = ggplot2::ggplot(data, ggplot2::aes_string(x = "time", y = y.name, group = ".algo", 
    linetype = ".algo", col = ".algo"))
  if (over.time == "dob")
    pl = pl + ggplot2::geom_point(size = 3)
  pl = pl + ggplot2::geom_line(data = mean.data, size = 1) 
  
  return(pl)
}

#' @title Plots Y traces of multiple optimization paths
#'
#' @description Plot function for \code{\link{renderYTraces}}
#'
#' @param opt.paths [\code{list}]\cr
#'   List of \code{OptPath} objects
#' @param over.time [\code{character}]\cr
#'   Should the traces be plotted versus the iteration number or the cumulated
#'   execution time? For the later, the opt.path has to contain a extra coloum
#'   names exec.time. Possible values are dob and exec.time, default is \code{dob}.
#'
#' @return [\code{NULL}]
#'
#' @export

plotYTraces = function(opt.paths, over.time = "dob") {
  print(renderYTraces(opt.paths, over.time))
  return(invisible(NULL))
}

