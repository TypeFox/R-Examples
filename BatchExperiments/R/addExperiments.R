#' @title Add experiemts to the registry.
#'
#' @description
#' Add experiments for running algorithms on problems
#' to the registry, so they can be executed later.
#'
#' @param reg [\code{\link{ExperimentRegistry}}]\cr
#'   Registry.
#' @param prob.designs [\code{character} | \code{\link{Design}} | list of \code{\link{Design}}]\cr
#'   Either problem ids, a single problem design or a list of problem designs,
#'   the latter two created by \code{\link{makeDesign}}.
#'   If missing, all problems are selected (without associating a design),
#'   and this is the default.
#' @param algo.designs [\code{character} | \code{\link{Design}} | list of \code{\link{Design}}]\cr
#'   Either algorithm ids, a single algorithm design or a list of algorithm designs,
#'   the latter two created by \code{\link{makeDesign}}.
#'   If missing, all algorithms are selected (without associating a design),
#'   and this is the default.
#' @param repls [\code{integer(1)}]\cr
#'   Number of replications.\cr
#'   Default is 1.
#' @param skip.defined [\code{logical}]\cr
#'   If set to \code{TRUE}, already defined experiments get skipped. Otherwise an error is thrown.\cr
#'   Default is FALSE.
#' @family add
#' @return Invisibly returns vector of ids of added experiments.
#' @examples
#' ### EXAMPLE 1 ###
#' reg = makeExperimentRegistry(id = "example1", file.dir = tempfile())
#'
#' # Define a problem:
#' # Subsampling from the iris dataset.
#' data(iris)
#' subsample = function(static, ratio) {
#'   n = nrow(static)
#'   train = sample(n, floor(n * ratio))
#'   test = setdiff(seq(n), train)
#'   list(test = test, train = train)
#' }
#' addProblem(reg, id = "iris", static = iris,
#'            dynamic = subsample, seed = 123)
#'
#' # Define algorithm "tree":
#' # Decision tree on the iris dataset, modeling Species.
#' tree.wrapper = function(static, dynamic, ...) {
#'   library(rpart)
#'   mod = rpart(Species ~ ., data = static[dynamic$train, ], ...)
#'   pred = predict(mod, newdata = static[dynamic$test, ], type = "class")
#'   table(static$Species[dynamic$test], pred)
#' }
#' addAlgorithm(reg, id = "tree", fun = tree.wrapper)
#'
#' # Define algorithm "forest":
#' # Random forest on the iris dataset, modeling Species.
#' forest.wrapper = function(static, dynamic, ...) {
#'   library(randomForest)
#'   mod = randomForest(Species ~ ., data = static, subset = dynamic$train, ...)
#'   pred = predict(mod, newdata = static[dynamic$test, ])
#'   table(static$Species[dynamic$test], pred)
#' }
#' addAlgorithm(reg, id = "forest", fun = forest.wrapper)
#'
#' # Define problem parameters:
#' pars = list(ratio = c(0.67, 0.9))
#' iris.design = makeDesign("iris", exhaustive = pars)
#'
#' # Define decision tree parameters:
#' pars = list(minsplit = c(10, 20), cp = c(0.01, 0.1))
#' tree.design = makeDesign("tree", exhaustive = pars)
#'
#' # Define random forest parameters:
#' pars = list(ntree = c(100, 500))
#' forest.design = makeDesign("forest", exhaustive = pars)
#'
#' # Add experiments to the registry:
#' # Use  previously defined experimental designs.
#' addExperiments(reg, prob.designs = iris.design,
#'                algo.designs = list(tree.design, forest.design),
#'                repls = 2) # usually you would set repls to 100 or more.
#'
#' # Optional: Short summary over problems and algorithms.
#' summarizeExperiments(reg)
#'
#' # Optional: Test one decision tree job and one expensive (ntree = 1000)
#' # random forest job. Use findExperiments to get the right job ids.
#' do.tests = FALSE
#' if (do.tests) {
#'   id1 = findExperiments(reg, algo.pattern = "tree")[1]
#'   id2 = findExperiments(reg, algo.pattern = "forest",
#'                          algo.pars = (ntree == 1000))[1]
#'   testJob(reg, id1)
#'   testJob(reg, id2)
#' }
#'
#' # Submit the jobs to the batch system
#' submitJobs(reg)
#'
#' # Calculate the misclassification rate for all (already done) jobs.
#' reduce = function(job, res) {
#'   n = sum(res)
#'   list(mcr = (n-sum(diag(res)))/n)
#' }
#' res = reduceResultsExperiments(reg, fun = reduce)
#' print(res)
#'
#' # Aggregate results using 'ddply' from package 'plyr':
#' # Calculate the mean over all replications of identical experiments
#' # (same problem, same algorithm and same parameters)
#' library(plyr)
#' vars = setdiff(names(res), c("repl", "mcr"))
#' aggr = ddply(res, vars, summarise, mean.mcr = mean(mcr))
#' print(aggr)
#'
#' \dontrun{
#' ### EXAMPLE 2 ###
#' # define two simple test functions
#' testfun1 = function(x) sum(x^2)
#' testfun2 = function(x) -exp(-sum(abs(x)))
#'
#' # Define ExperimentRegistry:
#' reg = makeExperimentRegistry("example02", seed = 123, file.dir = tempfile())
#'
#' # Add the testfunctions to the registry:
#' addProblem(reg, "testfun1", static = testfun1)
#' addProblem(reg, "testfun2", static = testfun2)
#'
#' # Use SimulatedAnnealing on the test functions:
#' addAlgorithm(reg, "sann", fun = function(static, dynamic) {
#'   upp = rep(10, 2)
#'   low = -upp
#'   start = sample(c(-10, 10), 2)
#'   res = optim(start, fn = static, lower = low, upper = upp, method = "SANN")
#'   res = res[c("par", "value", "counts", "convergence")]
#'   res$start = start
#'   return(res)
#' })
#'
#' # add experiments and submit
#' addExperiments(reg, repls = 10)
#' submitJobs(reg)
#'
#' # Gather informations from the experiments, in this case function value
#' # and whether the algorithm convergenced:
#' reduceResultsExperiments(reg, fun = function(job, res) res[c("value", "convergence")])
#' }
#' @aliases Experiment
#' @export
addExperiments = function(reg, prob.designs, algo.designs, repls = 1L, skip.defined = FALSE) {
  UseMethod("addExperiments")
}

#' @method addExperiments ExperimentRegistry
#' @export
addExperiments.ExperimentRegistry = function(reg, prob.designs, algo.designs, repls = 1L, skip.defined = FALSE) {
  checkExperimentRegistry(reg, strict = TRUE)
  BatchJobs:::syncRegistry(reg)

  # check prob.designs
  if (missing(prob.designs)) {
    prob.designs = lapply(dbGetAllProblemIds(reg), makeDesign)
  } else {
    if (is.character(prob.designs)) {
      prob.designs = lapply(prob.designs, makeDesign)
    } else if (is(prob.designs, "Design")) {
      prob.designs = list(prob.designs)
    } else if (is.list(prob.designs)) {
      checkListElementClass(prob.designs, "Design")
    } else {
      stop("Format of prob.designs not supported. Must be a character vector, a design or list of designs")
    }
    ids = unique(extractSubList(prob.designs, "id"))
    found = ids %in% dbGetAllProblemIds(reg)
    if (! all(found))
      stopf("%i problems have not been added to registry for designs: %s",
            sum(!found), collapse(ids[!found]))
  }

  # check algo.designs
  if (missing(algo.designs)) {
    algo.designs = lapply(dbGetAllAlgorithmIds(reg), makeDesign)
  } else {
    if (is.character(algo.designs)) {
      algo.designs = lapply(algo.designs, makeDesign)
    } else if (is(algo.designs, "Design")) {
      algo.designs = list(algo.designs)
    } else if (is.list(algo.designs)) {
      checkListElementClass(algo.designs, "Design")
    } else {
      stop("Format of algo.designs not supported. Must be a character vector, a design or list of designs")
    }
    ids = unique(extractSubList(algo.designs, "id"))
    found = ids %in% dbGetAllAlgorithmIds(reg)
    if (! all(found))
      stopf("%i algorithms have not been added to registry for designs: %s",
            sum(!found), collapse(ids[!found]))
  }

  repls = asCount(repls, positive = TRUE)
  assertFlag(skip.defined)

  f = function(xs) viapply(xs, function(x) x$designIter$n.states)
  n.exps = sum(outer(f(prob.designs), f(algo.designs)))
  info("Adding %i experiments / %i jobs to DB.", n.exps, n.exps*repls)
  if (n.exps == 0L)
    return(invisible(integer(0L)))

  # internal helper functions
  mq = function(lines, ..., con, bind.data = NULL) {
    q = sprintf(collapse(lines, sep = " "), ...)
    if(is.null(bind.data))
      return(dbGetQuery(con, q))
    return(dbGetPreparedQuery(con, q, bind.data = bind.data))
  }

  seripars = function(x) {
    rawToChar(serialize(x, connection = NULL, ascii = TRUE))
  }

  writeJobDefs = function(job.defs) {
    data = as.data.frame(do.call(rbind, lapply(job.defs, unlist)))
    mq("INSERT INTO tmp(prob_id, prob_pars, algo_id, algo_pars) VALUES(?, ?, ?, ?)",
       con = con, bind.data = data)
  }

  # establish persistent connection and create temporary table to fill
  # with job definitions
  con = BatchJobs:::dbConnectToJobsDB(reg, "rw")
  on.exit(dbDisconnect(con))

  # create temporary table for job definitions
  mq(c("CREATE TEMP TABLE tmp(job_def_id INTEGER, prob_id TEXT,",
       "prob_pars TEXT, algo_id TEXT, algo_pars TEXT)"), con = con)

  # write auxiliary temporary table with replication numbers
  mq("CREATE TEMP TABLE repls(repl INTEGER)", con = con)
  mq("INSERT INTO repls(repl) VALUES(?)",
     con = con, bind.data = data.frame(repl = seq_len(repls)))

  # create temporary view on cross product of repls and job_def_id
  mq(c("CREATE TEMP VIEW cp AS SELECT repls.repl, tmp.job_def_id FROM tmp",
       "CROSS JOIN repls"), con = con)



  # iterate to generate job definitions
  # write to temporary table every x definitions
  job.defs = BatchJobs:::buffer("list", 5000L, writeJobDefs)
  for (pd in prob.designs) {
    pd$designIter$reset()
    while (pd$designIter$hasNext()) {
      prob.pars = seripars(pd$designIter$nextElem())
      for (ad in algo.designs) {
        ad$designIter$reset()
        while (ad$designIter$hasNext()) {
          algo.pars = seripars(ad$designIter$nextElem())
          job.defs$push(list(prob_id = pd$id, prob_pars = prob.pars,
                             algo_id = ad$id, algo_pars = algo.pars))
        }
      }
    }
  }

  # add (remaining) defs to temporary job_defs table
  job.defs$clear()
  rm(job.defs)

  # query job_id to keep track of new ids
  max.job.id = mq("SELECT COALESCE(MAX(job_id), 0) AS x FROM %s_job_status", reg$id, con = con)$x

  # match for known job_def_id
  mq(c("UPDATE tmp SET job_def_id = (SELECT job_def_id FROM %s_job_def AS jd",
       "WHERE jd.prob_id = tmp.prob_id AND jd.algo_id = tmp.algo_id AND",
       "jd.prob_pars = tmp.prob_pars AND jd.algo_pars = tmp.algo_pars)"),
     reg$id, con = con)

  # test whether we would overwrite existing experiments
  if(!skip.defined) {
    if (mq("SELECT COUNT(job_def_id) AS n FROM tmp", con = con)$n > 0L)
      stop(paste("You have added identical experiments.",
                 "Either there are duplicated problem or algorithm ids or you have defined an experiment with the same parameters twice.",
                 "For the latter case use replications.",
                 "If you know what you're doing, look at skip.defined = TRUE.",
                 sep = "\n"))
  }

  # we start the transaction here, everything above is temporary
  dbBegin(con)
  ok = try({
    # insert new job defs
    mq(c("INSERT INTO %s_job_def(prob_id, prob_pars, algo_id, algo_pars)",
         "SELECT prob_id, prob_pars, algo_id, algo_pars FROM tmp",
         "WHERE job_def_id IS NULL"), reg$id, con = con)

    # update temporary table with new job defs
    mq(c("UPDATE tmp SET job_def_id = (SELECT job_def_id FROM %s_job_def AS jd WHERE",
         "jd.prob_id = tmp.prob_id AND jd.algo_id = tmp.algo_id AND",
         "jd.prob_pars = tmp.prob_pars AND jd.algo_pars = tmp.algo_pars)",
         "WHERE tmp.job_def_id IS NULL"), reg$id, con = con)

    # insert into job status table
    mq(c("INSERT INTO %1$s_job_status(job_def_id, repl)",
         "SELECT cp.job_def_id, cp.repl FROM cp WHERE NOT EXISTS",
         "(SELECT * FROM %1$s_job_status AS js WHERE",
         "cp.job_def_id = js.job_def_id AND cp.repl = js.repl)"),
       reg$id, con = con)


    # We could do this w/o bulk insert, but we are not allowed to
    # use external RNGs
    df = mq("SELECT job_id, pseed, repl FROM %s_expanded_jobs WHERE job_id > %i",
            reg$id, max.job.id, con = con)

    if(nrow(df) > 0L) {
      df$seed = BatchJobs:::addIntModulo(df$job_id, reg$seed - 1L)
      na = is.na(df$pseed)
      df$prob_seed[ na] = BatchJobs:::getRandomSeed(sum(na))
      df$prob_seed[!na] = BatchJobs:::addIntModulo(df$pseed[!na], df$repl[!na] - 1L)
      mq("UPDATE %s_job_status SET seed = ?, prob_seed = ? WHERE job_id = ?",
         reg$id, con = con, bind.data = df[c("seed", "prob_seed", "job_id")])
    }
  }, silent = TRUE)

  if(is.error(ok)) {
    dbRollback(con)
    errmsg = as.character(ok)
    # not really clean to match the english message here...
    # lets hope there are not localized versions of (R)SQLite out there
    if(grepl("UNIQUE constraint failed", errmsg, fixed = TRUE)) {
      stop(paste("You have added identical experiments.",
                 "Either there are duplicated problem or algorithm ids or you have defined an experiment with the same parameters twice.",
                 "For the latter case use replications.",
                 "If you know what you're doing, look at skip.defined = TRUE.",
                 sep = "\n"))
    } else {
      stopf("Error inserting new experiments: %s", errmsg)
    }
  }

  dbCommit(con)
  BatchJobs:::createShardedDirs(reg, df$job_id)
  invisible(df$job_id)
}
