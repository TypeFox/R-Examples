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

cmdLaunch <- function(verbose=FALSE) {
  if (is.na(verbose))
    verbose <- FALSE
  
  nwsName <- Sys.getenv('RSleighNwsName')
  userNwsName <- Sys.getenv('RSleighUserNwsName')
  nwsHost <- Sys.getenv('RSleighNwsHost')
  nwsPort <- as.integer(Sys.getenv('RSleighNwsPort'))
  maxWorkerCount <- as.integer(Sys.getenv('RSleighWorkerCount'))
  name <- Sys.getenv('RSleighName')

  launch(nwsName, nwsHost, nwsPort, maxWorkerCount, name, verbose, userNwsName)
}


launch <- function(nwsName, nwsHost, nwsPort, maxWorkerCount=-1,
      name=Sys.info()['nodename'], verbose=FALSE, userNwsName='__default') {
  nws <- netWorkSpace(nwsName, nwsHost, nwsPort, useUse=TRUE, create=FALSE)    
  userNws <- nwsUseWs(nws@server, userNwsName, create=FALSE)

  rank <- nwsFetch(nws, 'rankCount')
  if (rank < 0) {
    nwsStore(nws, 'rankCount', rank)
    stop('too late to join worker group')
  }
  else if ((maxWorkerCount >= 0) && (rank + 1 >= maxWorkerCount)) {
    nwsStore(nws, 'workerCount', maxWorkerCount)
    nwsStore(nws, 'rankCount', -1)
  }
  else {
    nwsStore(nws, 'rankCount', rank + 1)
  }

  # initialize for monitoring
  displayName = sprintf('%s@%d', name, rank)
  nwsDeclare(nws, displayName, 'single')
  nwsStore(nws, displayName, '0')
  nodeList = nwsFetch(nws, 'nodeList')
  if (nchar(nodeList) == 0)
    nodeList = displayName
  else
    nodeList = paste(nodeList, displayName)
  nwsStore(nws, 'nodeList', nodeList)

  # post some info about this worker
  logfile <- Sys.getenv('RSleighLogFile')
  names(logfile) <- NULL
  info <- Sys.info()
  nwsVersion <- paste(nwsPkgInfo(), collapse=' ')
  nwsStore(nws, 'worker info',
      list(host=info[['nodename']],
           os=info[['sysname']],
           pid=Sys.getpid(),
           R=R.version.string,
           nws=nwsVersion,
           rank=rank,
           logfile=logfile))

  # wait for all workers to join
  workerCount = nwsFind(nws, 'workerCount')

  # enter the main worker loop
  workerLoop(nws, displayName, rank, workerCount, verbose, userNws)

  # indicate exit.
  nwsStore(nws, 'bye', 1)
}
