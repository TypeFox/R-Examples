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

####
# sleighPending class
#
# represents a sleigh eachWorker/eachElem invocation in progress.
setClass('sleighPending',
         representation(nws='netWorkSpace', numTasks='numeric',
                        numSubmitted='numeric', accumulator='function',
                        barrierName='character', sleighState='environment',
                        state='environment'),
         prototype(nws=NULL))
setMethod('initialize', 'sleighPending',
function(.Object, nws, numTasks, numSubmitted, accumulator, bn, ss) {
  .Object@nws = nws
  .Object@numTasks = numTasks
  .Object@numSubmitted = numSubmitted
  .Object@accumulator = accumulator
  .Object@barrierName = bn
  .Object@sleighState = ss
  .Object@state = new.env()
  .Object@state$done = FALSE
  .Object
})

setMethod('show', 'sleighPending', function(object) {
  cat('\n')
  cat('NWS Sleigh Pending Object\n')
  show(object@nws)

  cat('Tasks submitted:', object@numSubmitted, '\n', sep='')

  status <- checkSleigh(object)
  if (status == 0)
    message <- 'Work completed.'
  else
    message <- paste(status, 'jobs still pending.')

  cat('Status:\t', message, '\n', sep='')
  cat('\n')
})

# return the number of results still outstanding.
setGeneric('checkSleigh', function(.Object) standardGeneric('checkSleigh'))
setMethod('checkSleigh', 'sleighPending',
function(.Object) {
  if (.Object@state$done) return (0) # could argue either way here... .

  v = nwsListVars(.Object@nws, showDataFrame=TRUE)
  n = tryCatch(v[v$Variable == 'result', 'NumValues'][[1]], error=function(e) 0)
  return (.Object@numSubmitted - n)
})

# collect all results.
#
# note: a lot of code is duplicated here and in the non-blocking sections of
# eachWorker and eachElem. refactor?
setGeneric('waitSleigh', function(.Object) standardGeneric('waitSleigh'))
setMethod('waitSleigh', 'sleighPending', function(.Object) {
  if (.Object@state$done) {
    stop('results already gathered.')
  }

  accum = if (is.null(body(.Object@accumulator))) NULL else .Object@accumulator
  val = if (is.null(accum)) vector('list', .Object@numTasks) else NULL
  accumargs = try(length(formals(accum)))  # results in 0 if accum is NULL

  if (.Object@numSubmitted > 0) {
    for (i in 1:.Object@numSubmitted) {
      repeat {
        r = nwsFetch(.Object@nws, 'result')
        # ignore everything but 'VALUE' messages
        if (is.list(r) && r$type == 'VALUE') break
      }

      # order results by rank for eachWorker, by tag for eachElem
      ind = if (.Object@barrierName != '') r$rank + 1 else r$tag

      if (is.null(accum)) {
        val[ind:(ind + length(r$value) - 1)] = r$value
      }
      else {
        if (accumargs == 0)
          try(accum())  # XXX should this be an error?
        else if (accumargs == 1)
          try(accum(r$value))
        else
          try(accum(r$value, ind:(ind + length(r$value) - 1)))
      }
    }
  }

  if (.Object@barrierName != '') {
    nwsStore(.Object@nws, .Object@barrierName, 1)
  }
  .Object@sleighState$occupied = FALSE
  .Object@sleighState$job = .Object@sleighState$job + 1
  .Object@state$done = TRUE

  val
})
