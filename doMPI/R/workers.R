#
# Copyright (c) 2009--2013, Stephen B. Weston
#
# This is free software; you can redistribute it and/or
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

#################
# master methods
#################

clusterSize <- function(cl, ...) {
  UseMethod('clusterSize')
}

closeCluster <- function(cl, ...) {
  UseMethod('closeCluster')
}

bcastSendToCluster <- function(cl, data, ...) {
  UseMethod('bcastSendToCluster')
}

sendToWorker <- function(cl, workerid, robj, ...) {
  UseMethod('sendToWorker')
}

recvFromAnyWorker <- function(cl, ...) {
  UseMethod('recvFromAnyWorker')
}

#################
# worker methods
#################

bcastRecvFromMaster <- function(cl, datalen, ...) {
  UseMethod('bcastRecvFromMaster')
}

sendToMaster <- function(cl, robj, ...) {
  UseMethod('sendToMaster')
}

recvFromMaster <- function(cl, ...) {
  UseMethod('recvFromMaster')
}

sinkWorkerOutput <- function(outfile) {
  if (outfile != "") {
    if (.Platform$OS.type == "windows" && outfile == "/dev/null")
      outfile <- "nul:"
    outcon <- file(outfile, open = "w")
    sink(outcon)
    sink(outcon, type = "message")
  }
}
