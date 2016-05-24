#
# Copyright (c) 2005-2008, REvolution Computing, Inc.
#
# NetWorkSpaces is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
# USA
#

# We use alternating barriers to synchronize eachWorker
# invocations. Their names are common to workers and sleighs.
barrierNames <- list('barrier0', 'barrier1')

# heuristic test for a closure
isClosure <- function(fun) {
  if (is.function(fun)) {
    e <- environment(fun)
    !is.null(e) && exists("environmentName", mode="function") &&
        identical(environmentName(e), "") &&
        length(ls(e, all.names=FALSE)) > 0
  } else {
    FALSE
  }
}

############################################################################
#  Sleigh code
#

sleigh <- function(...) {
  new("sleigh",...)
}

defaultSleighOptions <- new.env()

computeDefaultSleighOptions <- function(pkgpath) {
  # compute default value for scriptDir
  scriptDir = file.path(pkgpath, 'bin')

  # compute default value for nwsHost
  nwsHost = Sys.getenv('RSleighNwsHost')
  if (is.null(nwsHost) || (nchar(nwsHost) < 1))
    nwsHost = Sys.info()[['nodename']]

  # compute default value for nwsPort
  nwsPort = as.integer(Sys.getenv('RSleighNwsPort'))
  if (is.na(nwsPort))
    nwsPort = 8765

  # for some reason, R.home() uses backslashes on Windows,
  # even though .Platform says it should be forward slash.
  # this should fix the problem.
  rhome = gsub('\\\\', .Platform$file.sep, R.home())

  if (Sys.info()[['sysname']] == 'Windows') {
    scriptExec = scriptcmd
    scriptName = 'RNWSSleighWorker.py'
    workerWrapper = 'BackgroundLaunch.py'
    rprog = file.path(rhome, 'bin', 'Rterm.exe')
  }
  else {
    scriptExec = envcmd
    scriptName = 'RNWSSleighWorker.sh'
    workerWrapper = 'SleighWorkerWrapper.sh'  # XXX assumes ssh -f
    rprog = file.path(rhome, 'bin', 'R')
  }

  list(
      nwsHost = nwsHost,
      nwsHostRemote = NULL,
      nwsPort = nwsPort,
      nwsPortRemote = nwsPort,
      outfile = NULL,
      launch = 'local',
      workerCount = NULL,
      nodeList = NULL,
      scriptExec = scriptExec,
      wrapperDir = scriptDir, 
      scriptDir = scriptDir, 
      scriptName = scriptName,
      workerWrapper = workerWrapper, 
      workingDir = getwd(),
      logDir = NULL,
      user = NULL,
      passwd = NULL,
      simpleQuote = FALSE,
      wsNameTemplate = 'sleigh_ride_%04d',
      userWsNameTemplate = 'sleigh_user_%04d',
      verbose = FALSE,
      rprog = rprog,
      python = NULL
      )
}

####
# sleigh class
#
# represents a collection R processes running on a simple network of
# workstation pulling in tasks and generating results.


setClass('sleigh',
         representation(nodeList='character', nws='netWorkSpace',
                        userNws='netWorkSpace', nwsName='character',
                        userNwsName='character', nwss='nwsServer',
                        options='environment', state='environment'),
         prototype(nws=NULL, userNws=NULL, nwss=NULL))

setMethod('initialize', 'sleigh',
function(.Object, ...) {

  argList = list(...)

  # backward compatibility
  # check if nodeList is specified in the old way
  if (length(argList) > 0) {
    argName = names(argList[1])
    if (is.null(argName) || nchar(argName) == 0) {
      if (!is.vector(argList[[1]]))
        stop('argument 1 has no name and is not a vector')
      names(argList)[1] = 'nodeList'
      warning('nodeList should be passed using named variable, nodeList')
    }
  }

  # sanity check the optional arguments
  unrecog = names(argList)[!names(argList) %in% ls(defaultSleighOptions)]
  if (length(unrecog) > 0)
    stop('unused argument(s) ', paste(unrecog, collapse=', '))

  .Object@options = new.env()
  blendOptions(.Object@options, as.list(defaultSleighOptions))
  blendOptions(.Object@options, argList)
  opts = .Object@options

  .Object@state = new.env()
  .Object@state$bx = 1
  .Object@state$occupied = FALSE
  .Object@state$stopped = FALSE
  .Object@state$launch = opts$launch
  .Object@state$totalTasks = 0
  .Object@state$rankCount = 0
  .Object@state$job = 0

  if (!is.function(opts$launch) &&  !is.character(opts$launch)) {
    stop('unknown launch protocol.')
  }
  else if (is.character(opts$launch) &&
      !opts$launch %in% c('local', 'web', 'service')) {
    stop('unknown launch protocol.')
  }

  # set up the sleigh's netWorkSpace.
  .Object@nwss = new('nwsServer', serverHost=opts$nwsHost, port=opts$nwsPort)
  .Object@nwsName = nwsMktempWs(.Object@nwss, opts$wsNameTemplate)
  .Object@userNwsName = nwsMktempWs(.Object@nwss, opts$userWsNameTemplate)
  # NEED TO ADD ERROR HANDLING CODE HERE
  .Object@nws = nwsOpenWs(.Object@nwss, .Object@nwsName)
  .Object@userNws = nwsOpenWs(.Object@nwss, .Object@userNwsName)

  # initialize for monitoring
  nwsDeclare(.Object@nws, 'nodeList', 'single')
  nwsStore(.Object@nws, 'nodeList', '')
  nwsDeclare(.Object@nws, 'totalTasks', 'single')
  nwsStore(.Object@nws, 'totalTasks', '0')

  nwsDeclare(.Object@nws, 'rankCount', 'single')
  nwsStore(.Object@nws, 'rankCount', 0)
  nwsDeclare(.Object@nws, 'workerCount', 'single')
  nwsDeclare(.Object@nws, 'exported', 'fifo')

  if (is.function(opts$launch)) {
    if (is.null(opts$workerCount)) {
      .Object@nodeList = if (is.null(opts$nodeList))
                           rep('localhost', 3) else opts$nodeList
    }
    else {
      .Object@nodeList = if (is.null(opts$nodeList))
                           rep('localhost', opts$workerCount)
                         else
                           rep(opts$nodeList, length=opts$workerCount)
    }
    .Object@state$workerCount = length(.Object@nodeList)

    if (.Object@state$workerCount < 1) {
      close(.Object@nwss)
      stop('must have at least one worker in a sleigh')
    }

    for (i in 1:.Object@state$workerCount) {
      if (opts$verbose)
        opts$outfile = sprintf('%s_%04d.txt', .Object@nwsName, i)

      # since we are calling it in the constructor, maybe this cannot be
      # a method?
      addWorker(.Object@nodeList[i], .Object@nwsName, .Object@userNwsName,
                i, .Object@state$workerCount, opts)
    }
  }
  else if (opts$launch == 'service') {
    # remote launch using the "R Sleigh Service"
    service = tryCatch({
        nwsUseWs(.Object@nwss, 'RSleighService', create=FALSE)
      }, error=function(e) {
        close(.Object@nwss)
        stop('no sleigh services are running')
      })

    wsvars = nwsListVars(service, showDataFrame=TRUE)
    regworkers = unlist(wsvars$Variable, use.names=FALSE)
    user = if (is.null(opts$user)) Sys.info()[['login']] else opts$user

    # Note: we are only allowing execution on non-administrator sleigh services
    myworkers = regworkers[grep(paste('^', user, '@.', sep=''), regworkers)]

    if (!is.null(opts$nodeList)) {
      warning('ignoring user specified nodeList')
    }

    if (length(myworkers) < 1 || (!is.null(opts$workerCount) && opts$workerCount < 1)) {
      close(.Object@nwss)
      stop('must have at least one worker in a sleigh')
    }

    .Object@nodeList = if (!is.null(opts$workerCount))
                         rep(myworkers, opts$workerCount) else myworkers
    .Object@state$workerCount = length(.Object@nodeList)

    b = function(x) if (is.null(x)) '' else x
    for (i in 1:.Object@state$workerCount) {
      if (opts$verbose)
        opts$outfile = sprintf('%s_%04d.txt', .Object@nwsName, i)
      # XXX is '@' the best delimiter?
      request = sprintf('@%s@%s@%d@%d@%s@%s@%s@%s',
                        .Object@nwsName,
                        .Object@userNwsName,
                        .Object@state$workerCount,
                        i,
                        b(opts$workingDir),
                        b(opts$outfile),
                        b(opts$logDir),
                        user)
      if (opts$verbose)
        cat('command:', request, '\n')

      nwsStore(service, .Object@nodeList[i], request)
    }
  }
  else if (opts$launch == 'local') {
    # set workerCount if nobody else has
    if (is.null(opts$workerCount)) opts$workerCount <- 3
    if (opts$workerCount < 1) {
      close(.Object@nwss)
      stop('must have at least one worker in a sleigh')
    }

    .Object@state$workerCount = opts$workerCount
    .Object@nodeList = rep('localhost', opts$workerCount)
    for (i in 1:.Object@state$workerCount) {
      if (opts$verbose)
        opts$outfile = sprintf('%s_%04d.txt', .Object@nwsName, i)

      addWorker(.Object@nodeList[i], .Object@nwsName, .Object@userNwsName, i, .Object@state$workerCount, opts)
    }
  }
  else if (opts$launch == 'web') {
    cat(sprintf("Your ride is %s, don't forget 'DeleteMe...'.\n", .Object@nwsName))
    nwsStore(.Object@nws, 'runMe', sprintf("library(nws); launch('%s', '%s', %d, userNwsName='%s')",
                                           .Object@nwsName, opts$nwsHost, opts$nwsPort,
                                           .Object@userNwsName))

    tryCatch(nwsFetch(.Object@nws, 'deleteMeWhenAllWorkersStarted'), error=function(...) 0)
    nwsDeleteVar(.Object@nws, 'runMe')
    .Object@state$workerCount = nwsFetch(.Object@nws, 'rankCount')
    nwsStore(.Object@nws, 'workerCount', .Object@state$workerCount)
    nwsStore(.Object@nws, 'rankCount', -1)
    .Object@state$rankCount = -1
  }

  .Object
})

setMethod('show', 'sleigh', function(object) {
  cat('\n')
  cat('NWS Sleigh Object\n')
  show(object@nws)
  cat(object@state$workerCount, ' Worker Nodes:\t',
      paste(object@nodeList, collapse=', '), '\n', sep='')
  cat('\n')
})

# return a list that contains two values: workerCount and status
setGeneric('status', function(.Object, closeGroup=FALSE, timeout=0) standardGeneric('status'))
setMethod('status', 'sleigh',
function(.Object, closeGroup=FALSE, timeout=0) {
  if (.Object@state$rankCount < 0) {
    # join phase completed before
    return(list(numWorkers=workerCount(.Object), closed=1))
  }
  else {
    sleepTime = initialSleep = min(2.0, timeout)
    repeatedSleep = 1.0
    startTime = proc.time()[3]

    repeat {
      n = nwsFetch(.Object@nws, 'rankCount')
      if (n < 0) {
        # all workers joined
        if (.Object@options$verbose)
          cat('all ', workerCount(.Object), ' worker(s) have started\n')
        # this should not happened
        if (.Object@state$workerCount != nwsFind(.Object@nws, 'workerCount'))
          stop('Value for workerCount is not consistent.')
	.Object@state$rankCount = -1
        nwsStore(.Object@nws, 'rankCount', .Object@state$rankCount)
        return(list(numWorkers=workerCount(.Object), closed=1))
      }
      else {
        # three choices: close now, return not closed, or
        # sleep and try again
        if (proc.time()[3] - startTime >= timeout) {
          # time is up, so either close the join, or
          # return the current status
          if (closeGroup) {
            if (.Object@options$verbose)
              cat('closing group: ', n, ' workers\n')
            
            .Object@state$workerCount = n
            nwsStore(.Object@nws, 'workerCount', .Object@state$workerCount)
            .Object@state$rankCount = -1
            nwsStore(.Object@nws, 'rankCount', .Object@state$rankCount)
	    return(list(numWorkers=workerCount(.Object), closed=1))
          }
          else {
            if (.Object@options$verbose)
              cat('group not formed: ', n, ' worker(s)\n')

	    .Object@state$rankCount = n
            nwsStore(.Object@nws, 'rankCount', .Object@state$rankCount)
	    return(list(numWorkers=.Object@state$rankCount, closed=0))
          }
        }
        else {
          if (.Object@options$verbose)
            cat ('only ', n, ' worker(s): sleeping ...\n')

          nwsStore(.Object@nws, 'rankCount', n)
          Sys.sleep(sleepTime)
          sleepTime = repeatedSleep
        }
      }
    }
  }
})

if (! isGeneric('close'))
  setGeneric('close', function(con, ...) standardGeneric('close'))
setMethod('close', 'sleigh', function(con, ...) stopSleigh(con))

setGeneric('stopSleigh', function(.Object) standardGeneric('stopSleigh'))
setMethod('stopSleigh', 'sleigh', function(.Object) {
  if (.Object@state$stopped) return (invisible(NULL))

  if (!is.function(.Object@state$launch) &&
      identical(.Object@state$launch, 'web')) {
    nwsDeleteVar(.Object@nws, 'task')
  }
  else {
    nwsStore(.Object@nws, 'Sleigh ride over', 1)
  }
  Sys.sleep(3)
  exitCount = 0
  while (!is.null(nwsFetchTry(.Object@nws, 'bye'))) {
    exitCount = exitCount + 1
  }
  if (exitCount != .Object@state$workerCount) {
    cat(sprintf('Only %d of %d have exited.\n', exitCount, .Object@state$workerCount))
  }
  nwsDeleteWs(.Object@nwss, .Object@nwsName)
  close(.Object@nwss)
  .Object@state$stopped = TRUE
})


# run fun once on each worker of the sleigh. pass in a val from the
# range 1:#Workers
setGeneric('eachWorker',
           function(.Object, fun, ..., eo=NULL, DEBUG=FALSE) standardGeneric('eachWorker'))
setMethod('eachWorker', 'sleigh',
function(.Object, fun, ..., eo=NULL, DEBUG=FALSE) {
  if (DEBUG) browser()

  if (.Object@state$rankCount == -1 && .Object@state$workerCount < 1) {
    stop(paste('Worker group has been closed, and we have', .Object@state$workerCount, 'workers.'))
  }

  if (.Object@state$occupied) {
    stop('Sleigh is occupied')
  }

  if (.Object@state$stopped) {
    stop('Sleigh is stopped')
  }

  fun <- fun # need to force the argument (NJC: why?)

  nws = .Object@nws
  wc = .Object@state$workerCount

  blocking = TRUE
  accumulator = NULL
  closure = NULL
  if (!is.null(eo)) {
    if (is.environment(eo) || is.list(eo)) {
      if (!is.null(eo$blocking)) blocking = as.logical(eo$blocking)
      accumulator = eo$accumulator
      if (!is.null(eo$closure)) closure = as.logical(eo$closure)

      # check for unknown options
      if (is.list(eo)) {
        eo$blocking <- eo$accumulator <- eo$closure <- NULL
        if (length(eo) > 0)
          warning('ignoring unknown option(s): ',
            paste('"', names(eo), '"', sep='', collapse=', '))
      }
    }
    else {
      stop('options arg must be a list or environment.')
    }
  }

  # issue a warning if fun seems like a closure, and they aren't
  # explicitly enabled via the closure option.
  if (is.null(closure)) {
    closure <- TRUE
    if (isClosure(fun))
      warning('"fun" argument looks like a closure without enabling ',
        'via closure option', immediate.=TRUE)
  }

  # remove the enclosing environment of the function if closures are not
  # allowed.
  if (!closure)
    environment(fun) <- globalenv()

  # use alternating barrier to sync eachWorker invocations with the workers.
  bx = .Object@state$bx
  bn = barrierNames[[bx]]
  .Object@state$bx = bx%%2 + 1

  nwsFetchTry(.Object@nws, bn)

  # update the total number of submitted tasks
  .Object@state$totalTasks <- .Object@state$totalTasks + wc
  nwsStore(.Object@nws, 'totalTasks', as.character(.Object@state$totalTasks))

  # submit the tasks
  lapply(1:wc, storeTask, nws=nws, fun=fun, args=list(list(...)), barrier=TRUE,
         job=.Object@state$job)

  if (!blocking) {
    .Object@state$occupied = TRUE
    func = if (is.null(accumulator)) as.function(list(NULL)) else accumulator
    return (new('sleighPending', nws, wc, wc, func, bn, .Object@state))
  }

  val <- if (is.null(accumulator)) vector('list', wc) else NULL
  accumargs = try(length(formals(accumulator)))

  for (i in 1:wc) {
    repeat {
      r = nwsFetch(nws, 'result')
      # ignore everything but 'VALUE' messages
      if (is.list(r) && r$type == 'VALUE') break
    }
    if (is.null(accumulator)) {
      val[r$rank + 1] = r$value
    }
    else {
      if (accumargs == 0)
        accumulator()
      else if (accumargs == 1)
        accumulator(r$value)
      else
        accumulator(r$value, r$rank + 1)
    }
  }

  nwsStore(.Object@nws, bn, 1)

  # update the job id
  .Object@state$job = .Object@state$job + 1

  val
})


# run fun once for each element of a vector.
setGeneric('eachElem',
           function(.Object, fun, elementArgs=list(), fixedArgs=list(),
                    eo=NULL, DEBUG=FALSE) standardGeneric('eachElem'))
setMethod('eachElem', 'sleigh',
function(.Object, fun, elementArgs=list(), fixedArgs=list(), eo=NULL, DEBUG=FALSE) {
  if (DEBUG) browser()

  if (.Object@state$rankCount == -1 && .Object@state$workerCount < 1) {
    stop(paste('Worker group has been closed, and we have', .Object@state$workerCount, 'workers.'))
  }

  if (.Object@state$occupied) {
    stop('Sleigh is occupied')
  }

  if (.Object@state$stopped) {
    stop('Sleigh is stopped')
  }

  fun <- fun # need to force the argument (NJC: why?)

  nws = .Object@nws
  wc = .Object@state$workerCount

  argPermute = NULL
  blocking = TRUE
  lf = 0
  by = "row"
  chunkSize = 1
  accumulator = NULL
  elementFunc = NULL
  closure = NULL
  if (!is.null(eo)) {
    if (is.environment(eo) || is.list(eo)) {
      argPermute = eo$argPermute
      if (!is.null(eo$blocking)) blocking = as.logical(eo$blocking)
      if (!is.null(eo$loadFactor)) lf = as.numeric(eo$loadFactor)
      if (!is.null(eo$by)) by = match.arg(eo$by, c('row', 'column', 'cell'))
      if (!is.null(eo$chunkSize)) chunkSize = max(eo$chunkSize, 1)
      accumulator = eo$accumulator
      elementFunc = eo$elementFunc
      if (!is.null(eo$closure)) closure = as.logical(eo$closure)

      # check for unknown options
      if (is.list(eo)) {
        eo$argPermute <- eo$blocking <- eo$loadFactor <- eo$by <-
          eo$chunkSize <- eo$accumulator <- eo$elementFunc <- eo$closure <- NULL
        if (length(eo) > 0)
          warning('ignoring unknown option(s): ',
            paste('"', names(eo), '"', sep='', collapse=', '))
      }
    }
    else {
      stop('options arg must be a list or environment.')
    }
  }

  # issue a warning if fun seems like a closure, and they aren't
  # explicitly enabled via the closure option.
  if (is.null(closure)) {
    closure <- TRUE
    if (isClosure(fun))
      warning('"fun" argument looks like a closure without enabling ',
        'via closure option', immediate.=TRUE)
  }

  # remove the enclosing environment of the function if closures are not
  # allowed.
  if (!closure)
    environment(fun) <- globalenv()

  if (!is.list(elementArgs)) elementArgs = list(elementArgs)
  if (!is.list(fixedArgs)) fixedArgs = list(fixedArgs)

  if (length(elementArgs) > 0) {
    if (!is.null(elementFunc)) stop('elementFunc cannot be used with elementArgs')

    allTasks = unlist(lapply(elementArgs, countElement, by=by))
    # this allows for functions to be included, even though they aren't now
    numTasks = max(-1, allTasks, na.rm=TRUE)
    if (numTasks < 0) {
      numTasks = NA
    }
    else {
      # cat('got', numTasks, 'tasks\n')

      # check the length of the arguments
      for (taskLen in allTasks) {
        if (!is.na(taskLen) && numTasks %% taskLen != 0) {
          warning('elementArgs contains arguments of inconsistent length')
          break
        }
      }

      # update the total number of submitted tasks
      .Object@state$totalTasks <- .Object@state$totalTasks + numTasks
      nwsStore(.Object@nws, 'totalTasks', as.character(.Object@state$totalTasks))
    }
  }
  else if (!is.null(elementFunc)) {
    numTasks = NA
    nargs = length(formals(elementFunc))
    if (nargs > 2) stop('specified elementFunc function takes too many arguments')
    startingTasks = .Object@state$totalTasks
  }
  else {
    stop('either elementArgs or elementFunc must be specified')
  }

  if (blocking && lf > 0) {
    submitLimit = lf * wc
    if (submitLimit < wc) {
      submitLimit = wc
    }
  }
  else {
    submitLimit = Inf
  }

  allSubmitted = FALSE
  numSubmitted = 0
  numReturned = 0

  tag = 1
  currentTasks = 0

  val <- if (is.null(accumulator)) list() else NULL
  accumargs = try(length(formals(accumulator)))

  while (!allSubmitted || numReturned < numSubmitted) {
    if (!allSubmitted) {
      while (!allSubmitted && numSubmitted < submitLimit) {
        if (!is.null(elementFunc)) {
          argchunk = list()
          for (j in 1:chunkSize) {
            varArgs <- tryCatch(if (nargs == 0)
                elementFunc()
              else if (nargs == 1)
                elementFunc(currentTasks + 1)
              else
                elementFunc(currentTasks + 1, by), error=function(e) {
                  allSubmitted <<- TRUE
                  if (!is.null(e$message) && nchar(e$message) > 0)
                    warning(e$message)
                  NULL
                })
            if (allSubmitted) break
            if (!is.list(varArgs)) varArgs = list(varArgs)
            args = c(varArgs, fixedArgs)
            if (!is.null(argPermute)) args = args[argPermute]
            # cat('[function case] args:', paste(args, collapse=' '), '\n')
            argchunk[[j]] = args
            currentTasks = currentTasks + 1
          }
        }
        else {
          nTasks = min(numTasks - currentTasks, chunkSize)
          if (nTasks <= 0) {
            argchunk = list()
            allSubmitted = TRUE
          }
          else {
            if (nTasks > 1) {
              v = currentTasks:(currentTasks + nTasks - 1)
              varArgsChunk = lapply(1:length(elementArgs), function(j) {
                iv <- if (is.na(allTasks[j])) v + 1 else v %% allTasks[j] + 1
                getChunk(elementArgs[[j]], iv=iv, by=by)
              })
              argchunk = lapply(1:nTasks, function(i) {
                varArgs = lapply(varArgsChunk, getElement, i=i, by=by)
                args = c(varArgs, fixedArgs)
                if (!is.null(argPermute)) args = args[argPermute]
                # cat('[chunk case] args:', paste(args, collapse=' '), '\n')
                args
              })
            }
            else {
              varArgs = lapply(1:length(elementArgs), function(j) {
                i <- if(is.na(allTasks[j])) currentTasks + 1 else currentTasks %% allTasks[j] + 1
                getElement(elementArgs[[j]], i=i, by=by)
              })
              args = c(varArgs, fixedArgs)
              if (!is.null(argPermute)) args = args[argPermute]
              # cat('[standard case] args:', paste(args, collapse=' '), '\n')
              argchunk = list(args)
            }
            currentTasks = currentTasks + nTasks
          }
        }

        if (length(argchunk) > 0) {
          numSubmitted = numSubmitted + 1
          storeTask(nws, fun, argchunk, tag=tag, barrier=FALSE, job=.Object@state$job)
          tag = tag + length(argchunk)
        }
      }

      if (!is.null(elementFunc)) {
        # update the total number of submitted tasks
        .Object@state$totalTasks <- startingTasks + currentTasks
        nwsStore(.Object@nws, 'totalTasks', as.character(.Object@state$totalTasks))
      }

      if (!blocking) {
        .Object@state$occupied = TRUE
        # cat(sprintf('returning sleighPending object for %d tasks\n', nt))
        func = if (is.null(accumulator)) as.function(list(NULL)) else accumulator
        return (new('sleighPending', nws, currentTasks, numSubmitted, func, '', .Object@state))
      }
    }

    if (numReturned < numSubmitted) {
      repeat {
        r = nwsFetch(nws, 'result')
        # ignore everything but 'VALUE' messages
        if (is.list(r) && r$type == 'VALUE') break
      }
      if (is.null(accumulator)) {
        val[r$tag:(r$tag + length(r$value) - 1)] = r$value
      }
      else {
        if (accumargs == 0)
          accumulator()
        else if (accumargs == 1)
          accumulator(r$value)
        else
          accumulator(r$value, r$tag:(r$tag + length(r$value) - 1))
      }
      numReturned = numReturned + 1
      submitLimit = submitLimit + 1  # this can become > the number of tasks
    }
  }

  if (is.null(accumulator)) length(val) = currentTasks

  # update the job id
  .Object@state$job = .Object@state$job + 1

  val
})


setGeneric('rankCount', function(.Object) standardGeneric('rankCount'))
setMethod('rankCount', 'sleigh', function(.Object) .Object@state$rankCount)

setGeneric('workerCount', function(.Object) standardGeneric('workerCount'))
setMethod('workerCount', 'sleigh', function(.Object) .Object@state$workerCount)

setGeneric('netWorkSpaceObject', function(.Object) standardGeneric('netWorkSpaceObject'))
setMethod('netWorkSpaceObject', 'sleigh', function(.Object) .Object@nws)

wsVarName <- function(name, worker) {
  if (is.null(worker)) {
    sprintf('env_%s', name)
  } else {
    sprintf('env_%d_%s', worker, name)
  }
}

setGeneric('export',
           function(.Object, xName, xVal, worker=NULL) standardGeneric('export'))
setMethod('export', 'sleigh',
function(.Object, xName, xVal, worker=NULL) {
  # sleigh error checking
  if (.Object@state$occupied) stop('Sleigh is occupied')
  if (.Object@state$stopped) stop('Sleigh is stopped')

  # argument error checking
  if (missing(xName)) stop('no value specified for xName argument')
  if (missing(xVal)) stop('no value specified for xVal argument')
  if (! is.character(xName)) stop('xName must be a character variable')
  if (! is.null(worker)) {
    if (! is.numeric(worker)) stop('worker value must be numeric')
    if (length(worker) > 1) stop('only one worker can be specified')
    if (worker < 0) stop('worker value must be positive')
    if (worker >= .Object@state$workerCount)
      stop('worker value is too large for this sleigh')
  }

  wsVar <- wsVarName(xName, worker)
  nwsDeclare(.Object@nws, wsVar, 'single')
  nwsStore(.Object@nws, wsVar, xVal)
  nwsStore(.Object@nws, 'exported', list(worker=worker, name=xName, wsVar=wsVar))
  invisible(NULL)
})

setGeneric('unexport',
           function(.Object, xName, worker=NULL) standardGeneric('unexport'))
setMethod('unexport', 'sleigh',
function(.Object, xName, worker=NULL) {
  # sleigh error checking
  if (.Object@state$occupied) stop('Sleigh is occupied')
  if (.Object@state$stopped) stop('Sleigh is stopped')

  # argument error checking
  if (missing(xName)) stop('no value specified for xName argument')
  if (! is.character(xName)) stop('xName must be a character variable')
  if (! is.null(worker)) {
    if (! is.numeric(worker)) stop('worker value must be numeric')
    if (length(worker) > 1) stop('only one worker can be specified')
    if (worker < 0) stop('worker value must be positive')
    if (worker >= .Object@state$workerCount) stop('worker value is too large')
  }

  wsVar <- wsVarName(xName, worker)
  tryCatch({
    nwsDeleteVar(.Object@nws, wsVar)
  }, error=function(e) {
    stop('cannot unexport a variable that was not exported: ', wsVar)
  })
  nwsStore(.Object@nws, 'exported', list(worker=worker, name=xName, wsVar=NULL))
  invisible(NULL)
})

setGeneric('workerInfo',
           function(.Object) standardGeneric('workerInfo'))
setMethod('workerInfo', 'sleigh',
function(.Object) {
  n <- .Object@state$workerCount
  host <- as.character(rep(NA, n))
  os <- as.character(rep(NA, n))
  pid <- as.integer(rep(NA, n))
  R <- as.character(rep(NA, n))
  nws <- as.character(rep(NA, n))
  rank <- as.integer(rep(NA, n))
  logfile <- as.character(rep(NA, n))

  it <- nwsIFindTry(.Object@nws, 'worker info')
  x <- it()
  while (! is.null(x)) {
    i <- as.integer(x$rank) + 1
    host[i] <- x$host
    os[i] <- x$os
    pid[i] <- as.integer(x$pid)
    R[i] <- x$R
    nws[i] <- x$nws
    rank[i] <- as.integer(x$rank)
    logfile[i] <- x$logfile
    x <- it()
  }
  data.frame(host=host, os=os, pid=pid, R=R, nws=nws,
             rank=rank, logfile=logfile, stringsAsFactors=FALSE)
})
