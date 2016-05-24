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

startnode <- function(label, parent=NULL) {
  eparent <- if (is.null(parent)) emptyenv() else parent
  node <- new.env(parent=eparent)  # XXX not sure if there is any point
  node$i <- 1  # initialize my own "child index"
  node$label <- label
  node$starttime <- proc.time()[[3]]

  if (!is.null(parent)) {
    if (!exists('starttime', where=parent, inherits=FALSE))
      warning(sprintf('parent node %s has not been started', parent$label))
    if (exists('finishtime', where=parent, inherits=FALSE))
      warning(sprintf('parent node %s has already finished', parent$label))
    assign(getchildname(parent$i), node, pos=parent)
    parent$i <- parent$i + 1
  }

  node
}

finishnode <- function(node, newlabel=NULL) {
  if (!exists('starttime', where=node, inherits=FALSE))
    stop(sprintf('node %s has not been started', node$label))
  if (exists('finishtime', where=node, inherits=FALSE))
    stop(sprintf('node %s has already finished', node$label))
  node$finishtime <- proc.time()[[3]]
  if (!is.null(newlabel))
    node$label <- newlabel
  invisible(node)
}

getchildname <- function(i) {
  paste('node', i, sep='_')
}

getchild <- function(node, i) {
  get(getchildname(i), pos=node, inherits=FALSE)
}

nodetime <- function(node) {
  if (!exists('starttime', where=node, inherits=FALSE))
    stop('nodetime called on an invalid node')
  if (!exists('finishtime', where=node, inherits=FALSE))
    stop(sprintf('node %s has not finished', node$label))
  node$finishtime - node$starttime
}

displaynode <- function(node, recurse=TRUE, level=0) {
  total <- nodetime(node)
  childtime <- 0
  numchildren <- node$i - 1

  for (i in seq(length=numchildren)) {
    childtime <- childtime + nodetime(getchild(node, i))
  }

  indent <- paste(rep('    ', level), collapse='')
  cat(sprintf('%s%s: total time: %f; self time: %f\n',
              indent, node$label, total, total - childtime))

  if (recurse) {
    for (i in seq(length=numchildren)) {
      displaynode(getchild(node, i), recurse=recurse, level=level+1)
    }
  }
  invisible(NULL)
}
