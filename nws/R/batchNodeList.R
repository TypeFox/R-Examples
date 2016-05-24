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

batchNodeList <- function() {
  if (Sys.getenv('PE_HOSTFILE') != '') {
    sgeNodeList()
  }
  else if (Sys.getenv('LSB_HOSTS') != '') {
    lsfNodeList()
  }
  else if (Sys.getenv('PBS_NODEFILE') != '') {
    pbsNodeList()
  }
  else {
    stop('cannot determine the kind of batch queueing system used')
  }
}

sgeNodeList <- function() {
  hostfile <- Sys.getenv('PE_HOSTFILE')
  if (hostfile == '')
    stop('environment variable PE_HOSTFILE is not defined')

  nodeList <- readLines(hostfile)
  nodeList <- strsplit(nodeList, '[[:space:]]+')

  tryCatch({
      unlist(lapply(nodeList, function(x) rep(x[1], as.integer(x[2]))))
    }, error=function(e) {
      stop('hostfile has bad format: ', hostfile)
    })
}

lsfNodeList <- function() {
  hostlist <- Sys.getenv('LSB_HOSTS')
  if (hostlist == '')
    stop('environment variable LSB_HOSTS is not defined')

  nodeList <- unlist(strsplit(hostlist, '[[:space:]]+'), use.names=FALSE)
  x <- grep('^[[:space:]]*$', nodeList)
  if (length(x) > 0)
    nodeList <- nodeList[-x]
  nodeList
}

pbsNodeList <- function() {
  hostfile <- Sys.getenv('PBS_NODEFILE')
  if (hostfile == '')
    stop('environment variable PBS_NODEFILE is not defined')

  readLines(hostfile)
}
