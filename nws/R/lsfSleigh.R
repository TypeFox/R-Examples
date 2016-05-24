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

# This function will properly check to see if the R script has been
# launched in a parallel fashion.  If so, it will start jobs on the
# provided nodes using ssh.  If not, it expects that the parameter
# 'N' contains the number of nodes to start via 'bsub'.  Extra
# arguments to bsub can be provided using the lsfoptions parameter.

lsfSleigh <- function(n, lsfOptions=c(), ...) {
  hostList <- unlist(Sys.getenv("LSB_HOSTS"))
  if(! is.null(hostList) && hostList > " ") {
    hosts <- unlist(strsplit(hostList, " "))
    hosts <- hosts[hosts > " "]
    if(! missing(n) && n != length(hosts)) {
      warning("Number of requested nodes (", n,
              ") does not match number of nodes assigned by LSF (",
              length(hosts),
              ").  Using the ", length(hosts), " hosts LSF has assigned.")
    }
    n <- length(hosts)
    cmd <- sshcmd
  }
  else if(missing(n)) {
    stop("No requested and no hosts assigned by LSF")
  }
  else {
    hosts <- as.character(1:n)
    cmd <- function(...) {
      c('bsub', lsfOptions)
    }
  }

  sleigh(launch=cmd, nodeList=hosts, ...)
}
