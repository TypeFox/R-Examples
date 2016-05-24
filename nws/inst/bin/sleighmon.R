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

library(nws)
library(lattice)

quote_plus <- function(x) {
  x <- gsub(';', '%3B', x, fixed=TRUE)
  x <- gsub('/', '%2F', x, fixed=TRUE)
  x <- gsub('?', '%3F', x, fixed=TRUE)
  x <- gsub(':', '%3A', x, fixed=TRUE)
  x <- gsub('@', '%40', x, fixed=TRUE)
  x <- gsub('&', '%26', x, fixed=TRUE)
  x <- gsub('=', '%3D', x, fixed=TRUE)
  x <- gsub('+', '%2B', x, fixed=TRUE)
  x <- gsub('$', '%24', x, fixed=TRUE)
  x <- gsub(',', '%2C', x, fixed=TRUE)
  x <- gsub(' ', '+', x, fixed=TRUE)
  x
}

escape <- function(x) {
  x <- gsub('&', '&amp;', x, fixed=TRUE)
  x <- gsub('<', '&lt;', x, fixed=TRUE)
  x <- gsub('>', '&gt;', x, fixed=TRUE)
  x
}

monitor <- function(ws, filename) {
  nodelist <- nwsFind(ws, 'nodeList')
  nodes <- unlist(strsplit(nodelist, " "))
  tasks <- unlist(lapply(nodes, function(n) as.integer(nwsFind(ws, n))))
  labels <- sub('localhost@', '', nodes)
  sumTasks <- sum(tasks)
  totalTasks <- max(as.integer(nwsFind(ws, 'totalTasks')), sumTasks)
  xlim <- c(0, max(totalTasks, 1))
  main <- nwsFindTry(ws, 'MainTitle', 'Sleigh Monitor')
  sub <- paste('Completed', sumTasks, 'of', totalTasks, 'Tasks')
  sub <- nwsFindTry(ws, 'SubTitle', sub)

  h <- 8  # XXX ?
  res <- 72
  pixels <- res * h

  tryCatch(png(filename=filename, height=pixels, width=576, pointsize=12),
           error=function(e) bitmap(file=filename, height=h, width=8, res=res, pointsize=12))

  x <- barchart(labels~tasks,
                horizontal=TRUE,
                xlab='Tasks Executed',
                ylab='Workers',
                xlim=xlim,
                col=rainbow(length(labels)),
                main=main)

  total <- 'Total'
  y <- barchart(total~sumTasks,
                horizontal=TRUE,
                xlab=NULL,
                ylab=NULL,
                xlim=xlim,
                scales=list(draw=FALSE),
                col='red',
                sub=sub)

  print(x, position=c(0,0.1,1,1), more=TRUE)
  print(y, position=c(0,0,1,0.15))

  dev.off()
}

getErrLog <- function(vws) {
  refresh <- c('doit?op=showMonitor&monName=Sleigh+Monitor&wsName=',
               quote_plus(nwsWsName(vws)))
  refresh <- paste(refresh, collapse='')
  it <- nwsIFindTry(vws, 'logError')
  x <- list()
  tmp <- '        <tr class="%s"><td>%s</td></tr>'
  oe <- c('even', 'odd')
  j <- 0
  while (!is.null(i <- it())) {
    x <- c(x, list(sprintf(tmp, oe[j %% 2 + 1], escape(i))))
    j <- j + 1
  }

  if (length(x) > 0) {
    hdr1 <- '<html>
<head>
  <title>Sleigh Worker Error Log</title>
  <style>
    ul.menu {
      display: block;
      padding: 0 0.3em 0.3em 0.3em;
      background-color: #777;
      margin: 0 0 1em 0;
    }
    ul.menu li {
      display: inline;
      padding: 0 0 0 1em;
    }
    ul.menu a:hover {
      background-color: black;
      text-decoration: none;
    }
    ul.menu a:link, ul.menu a:visited {
      color: white;
      text-decoration: none;
    }
    a:link, a:visited {
      color: blue;
      text-decoration: none;
    }
    .tableheader {
      background-color: #cecece;
    }
    .even {
      background-color: #eee;
    }
    .odd {
      background-color: #dedede;
    }
    .error {
      color: #EE1111;
      font-weight: bold;
      margin: 20;
    }
    .confirm {
      color: #EE1111;
      font-weight: bold;
      margin: 20;
    }
    .info {
      color: black;
      margin: 20;
    }
    body {
      border: 0;
      padding: 0;
      margin: 0;
      background-color: #efefef;
      font-family: helvetica, sans-serif;
    }
    h1 {
      padding: 0.3em;
      background-color: #777;
      color: white;
      font-size: 1.6em;
      font-style: italic;
      font-weight: bold;
      margin-bottom: 0;
    }
    th {
      text-align: left;
    }
    table form {
      margin-bottom: 0;
    }
    input[value=X] {
      color: #EE1111;
      font-weight: bold;
    }
    .nwstable {
      margin: 0 0 0 1em;
    }

    input {
      display: block;
      width: 50em;
      float: left;
      margin-bottom: 10px;
      margin-top: 1em;
      margin-left: 1.5em;
    }
    label {
      display: block;
      text-align: right;
      float: left;
      width: 75px;
      padding-right: 20px;
      margin-top: 1em;
    }
    br {
      clear: left;
    }
    .buttonSubmit {
      width: 75px;
      margin-left: 15px;
    }
    .hidden {
      display: none;
    }
  </style>
  <body>
    <h1>Sleigh Worker Error Log</h1>
    <ul class="menu">
      <li><a href="doit?op=listWss"> NetWorkSpaces </a></li>
      <li><a href="'

    hdr2 <- '"> Refresh </a></li>
    </ul>
    <table cellpadding="4" class="nwstable">
      <tbody>'

    ftr <- '      </tbody>
    </table>
  </body>
</html>
'
    c(paste(hdr1, refresh, hdr2, sep=''), x, ftr)
  } else {
    x
  }
}

host <- 'localhost'
port <- 8765

gotargs <- FALSE
args <- commandArgs()

while (length(args) > 0) {
  a <- args[1]
  args <- args[-1]

  if (!gotargs) {
    if (a == '--args') gotargs <- TRUE
  } else {
    a <- match.arg(a, c('-host', '-port'))
    if (length(args) == 0) stop('option ', a, ' takes a required argument')
    v <- args[1]
    args <- args[-1]
    assign(substring(a, 2), switch(a, '-host' = v, '-port' = as.integer(v)))
  }
}

bws <- netWorkSpace('Sleigh Monitor', host, port)

reqnum <- 0
while (TRUE) {
  # cat('waiting for request\n')
  vwsname <- NA
  replyVarName <- 'reply'
  a <- nwsFetch(bws, 'request')
  # cat('number of arguments', a, '\n')
  numargs <- as.integer(a)
  if (!is.na(numargs) && numargs > 0) {
    for (i in 1:numargs) {
      # cat('waiting for argument', i, '\n')
      val <- nwsFetch(bws, 'request')
      # cat('argument:', val, '\n')
      v <- strsplit(val, "")
      x <- grep("=", v[[1]])
      if (length(x) > 0 && x[1] > 1) {
        name <- substring(val, 1, x[1] - 1)
        # cat('arg name:', name, '\n')
        value <- substring(val, x[1] + 1)
        # cat('arg value:', value, '\n')
        if (name == 'wsName') {
          vwsname <- value
        } else if (name == 'replyVarName') {
          replyVarName <- value
        }
      } else {
        cat('bad argument:', val)
        vwsname <- NA
        break
      }
    }
  }

  if (!is.na(vwsname)) {
    # cat('got request for workspace', vwsname, '\n')

    # open the workspace and create an image file from it
    vws <- nwsUseWs(bws@server, vwsname)  # XXX should use 'create=FALSE'
    # cat('opening sleigh workspace\n')

    if (length(errs <- getErrLog(vws)) == 0) {
      reqnum <- reqnum + 1
      imagefile <- sprintf('mon_%d.png', reqnum)
      monitor(vws, imagefile)
      # cat('returned from monitor function\n')

      # read the image file into a variable
      size <- file.info(imagefile)$size
      # cat(imagefile, 'has', size, 'bytes\n')
      image <- readBin(imagefile, 'raw', size)

      # send the image data to the web interface and remove the image file
      # cat('sending reply', length(image), 'bytes long\n')
      nwsStore(bws, replyVarName, '3')
      nwsStore(bws, replyVarName, 'content-type=image/png')
      nwsStore(bws, replyVarName, 'cache-control=no-cache')
      nwsStore(bws, replyVarName, 'refresh=5')
      nwsStore(bws, replyVarName, image)
      # cat('removing image file\n')
      file.remove(imagefile)
    } else {
      # display error message that have been logged to the workspace
      nwsStore(bws, replyVarName, '1')
      nwsStore(bws, replyVarName, 'cache-control=no-cache')
      nwsStore(bws, replyVarName, do.call(paste, c(errs, list(''), list(sep='\n'))))
    }
  } else {
    cat('bad request\n')
    nwsStore(bws, replyVarName, '2')
    nwsStore(bws, replyVarName, 'content-type=text/plain')
    nwsStore(bws, replyVarName, 'cache-control=no-cache')
    nwsStore(bws, replyVarName, 'wsName not specified')
  }
}
