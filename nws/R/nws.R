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

# this is a manifest constant, how to approriately handle?
nwsRFP = 3*2^24

# a utility function to read exactly n bytes from the socket connection.
nwsRecvN <- function(s, n, rawflag=FALSE) {
  if (n > 0) {
    b = readBin(s, what='raw', n=n)
    m = length(b)
    n = n - m
    if (n <= 0) return(if (rawflag) b else rawToChar(b))
    if (m == 0) stop("failed to read from nws socket")

    # we didn't get all of the data, so save the raw vector in a list
    # that we'll concatenate when we do have it all
    rlen = 50
    r = vector('list', rlen)
    i = 1
    r[i] = list(b)

    repeat {
      b = readBin(s, what='raw', n=n)
      i = i + 1
      if (i > rlen) {
        # we ran out of space in our list, so double its length
        rlen = 2 * rlen
        length(r) = rlen
      }
      r[i] = list(b)
      m = length(b)
      n = n - m
      if (n <= 0) break
      if (m == 0) stop("failed to read from nws socket")
    }

    # truncate the list, concatenate the raw vectors,
    # and convert to a single character string
    length(r) = i
    return(if (rawflag) do.call(c, r) else rawToChar(do.call(c, r)))
  }
  else {
    return(if (rawflag) raw(0) else '')
  }
}

nwsServer <- function(...) {
  new("nwsServer", ...)
}

# class respresenting connection to a netWorkSpace server.
setClass('nwsServer', representation(nwsSocket='ANY', port='numeric',
         serverHost='character', cookieProtocol='logical'))

setMethod('initialize', 'nwsServer',
          function(.Object, serverHost='localhost', port=8765) {
            .Object@serverHost = serverHost
            .Object@port = port

            if (Sys.info()[['sysname']] == 'Windows') {
              # on windows, socketConnection will wait for the full timeout,
              # even if no one is listening on the specified server port.
              # make.socket doesn't, so we'll use it to throw an exception
              # if no one is listening.
              tmpsock = make.socket(serverHost, port)
              close.socket(tmpsock)
            }

            # temporarily change the timeout while creating the socketConnection.
            # we will block for up to a year for data on this socket.
            old.timeout = options(timeout=32140800)[[1]]
            .Object@nwsSocket = tryCatch(socketConnection(serverHost, port=port, open='a+b', blocking=TRUE),
                                         finally=options(timeout=old.timeout))

            # tell the server that we're a new client
            writeBin(charToRaw('1112'), .Object@nwsSocket)
            handshake <- nwsRecvN(.Object@nwsSocket, 4)
            .Object@cookieProtocol <- handshake != '2222'
            .Object
          })

setGeneric('nwsDeleteWs', function(.Object, wsName) standardGeneric('nwsDeleteWs'))
setGeneric('nwsListWss', function(.Object, showDataFrame=FALSE) standardGeneric('nwsListWss'))
setGeneric('nwsMktempWs', function(.Object, wsNameTemplate='__Rws__%010d') standardGeneric('nwsMktempWs'))
setGeneric('nwsOpenWs', function(.Object, wsName, space=NULL, ...) standardGeneric('nwsOpenWs'))
setGeneric('nwsUseWs', function(.Object, wsName, space=NULL, ...) standardGeneric('nwsUseWs'))

setMethod('nwsDeleteWs', 'nwsServer',
          function(.Object, wsName) {
            op = 'delete ws'
            s = .Object@nwsSocket

            writeBin(charToRaw(sprintf('0002%020d%s%020d%s', nchar(op), op, nchar(wsName), wsName)), s)

            # status, unused at the moment.
            bb = nwsRecvN(s, 4)
          })

setMethod('nwsListWss', 'nwsServer',
          function(.Object, showDataFrame=FALSE) {
            op = 'list wss'
            s = .Object@nwsSocket
            writeBin(charToRaw(sprintf('0001%020d%s', nchar(op), op)), s)

            status = as.integer(nwsRecvN(s, 4))
            desc = nwsRecvN(s, 20)
            if (.Object@cookieProtocol)
              cookie <- nwsRecvN(s, 40)

            ret <- nwsRecvN(s, as.integer(nwsRecvN(s, 20)))
            if (showDataFrame==FALSE)
              ret
            else {
              ## convert response into an R data frame
              ret <- unlist(strsplit(ret, "\n"))
              retval <- list()
              fields <- list()
              i = 1
              while (i <= length(ret)) {
                line <- unlist(strsplit(ret[i], "\t"))

                # convert each field to correct type
                fields[1] = FALSE
                if (substr(line[1], 1, 1)=='>')
                  fields[1] = TRUE
                fields[2] = substr(line[1], 2, nchar(line[1]))  # workspace name
                fields[3] = line[2]
                fields[4] = as.logical(line[3])
                fields[5] = as.integer(line[4])
                if (is.na(line[5]))
                  fields[6] = ""
                else
                  fields[6] = line[5]

                retval = c(retval, list(fields))
                i = i+1
              }

              if (length(retval) > 0) {
                names(retval) <- seq(along=retval)
                retval <- do.call(rbind, retval)
                colnames(retval) <-
                  c("Owned", "Name", "Owner", "Persistent", "NumVariables", "Variables")
              }
              retval <- data.frame(retval)
              retval
            }

          })

setMethod('nwsMktempWs', 'nwsServer',
          function(.Object, wsNameTemplate) {
            if (!is.character(wsNameTemplate))
              stop('workspace name must be a string')

            op = 'mktemp ws'
            s = .Object@nwsSocket
            writeBin(charToRaw(sprintf('0002%020d%s%020d%s',
                nchar(op), op, nchar(wsNameTemplate), wsNameTemplate)), s)
            status = as.integer(nwsRecvN(s, 4))
            desc = nwsRecvN(s, 20) # unused at the moment.
            if (.Object@cookieProtocol)
              cookie <- nwsRecvN(s, 40) # unused at the moment.
            n <- as.integer(nwsRecvN(s, 20))
            name <- nwsRecvN(s, n)
            if (status != 0) stop('mktempWs failed')
            name
          })

setMethod('nwsOpenWs', 'nwsServer',
          function(.Object, wsName, space=NULL, ...) {
            # sanity check the optional arguments
            opts = list(...)
            unrecog = names(opts)[!names(opts) %in% c('create', 'persistent')]
            if (length(unrecog) > 0)
              stop('unused argument(s) ', paste(unrecog, collapse=', '))

            # if invoked directly by user, we need to create a space
            # instance. if invoked via networkspace constructor, use the
            # space passed in.
            if (is.null(space)) {
              serverWrap = new.env()
              serverWrap$server = .Object
              space = new('netWorkSpace', wsName=wsName, serverWrap=serverWrap)
            }

            op = 'open ws'
            owner = sprintf('%d', Sys.getpid())
            p = 'no'
            if (!is.null(opts$persistent) && opts$persistent) p = 'yes'

            s = .Object@nwsSocket

            if (is.null(opts$create) || opts$create) {
              writeBin(charToRaw(sprintf('0004%020d%s%020d%s%020d%s%020d%s',
                                      nchar(op), op,
                                      nchar(wsName), wsName,
                                      nchar(owner), owner,
                                      nchar(p), p)), s)
            }
            else {
              create = 'no'
              writeBin(charToRaw(sprintf('0005%020d%s%020d%s%020d%s%020d%s%020d%s',
                                      nchar(op), op,
                                      nchar(wsName), wsName,
                                      nchar(owner), owner,
                                      nchar(p), p,
                                      nchar(create), create)), s)
            }

            status = as.integer(nwsRecvN(s, 4))
            if (status != 0) stop(paste("workspace", wsName, "doesn't exist"))
            space
          })

setMethod('nwsUseWs', 'nwsServer',
          function(.Object, wsName, space=NULL, ...) {
            # sanity check the optional arguments
            opts = list(...)
            unrecog = names(opts)[!names(opts) %in% c('create', 'persistent')]
            if (length(unrecog) > 0)
              stop('unused argument(s) ', paste(unrecog, collapse=', '))

            # see nwsOpenWs
            if (is.null(space)) {
              serverWrap = new.env()
              serverWrap$server = .Object
              space = new('netWorkSpace', wsName=wsName, serverWrap=serverWrap)
            }

            op = 'use ws'
            owner = ''
            p = 'no'

            s = .Object@nwsSocket

            if (is.null(opts$create) || opts$create) {
              writeBin(charToRaw(sprintf('0004%020d%s%020d%s%020d%s%020d%s',
                                      nchar(op), op,
                                      nchar(wsName), wsName,
                                      nchar(owner), owner,
                                      nchar(p), p)), s)
            }
            else {
              create = 'no'
              writeBin(charToRaw(sprintf('0005%020d%s%020d%s%020d%s%020d%s%020d%s',
                                      nchar(op), op,
                                      nchar(wsName), wsName,
                                      nchar(owner), owner,
                                      nchar(p), p,
                                      nchar(create), create)), s)
            }

            status = as.integer(nwsRecvN(s, 4))
            if (status != 0) stop(paste("workspace", wsName, "doesn't exist"))
            space
          })

if (!isGeneric('close'))
  setGeneric('close', function(con, ...) standardGeneric('close'))
setMethod('close', 'nwsServer', function(con, ...) close(con@nwsSocket))

netWorkSpace <- function(...) {
  new("netWorkSpace", ...)
}

# class representing a netWorkSpace.
setClass('netWorkSpace', representation(server='nwsServer',
         wsName='character', cookieProtocol='logical'),
         prototype(server=NULL))

setMethod('initialize', 'netWorkSpace',
          function(.Object, wsName='__default', serverHost='localhost',
                   port=8765, useUse=FALSE, serverWrap=NULL, ...) {
            # sanity check the optional arguments
            argNames = names(list(...))
            unrecog = argNames[!argNames %in% c('create', 'persistent')]
            if (length(unrecog) > 0)
              stop('unused argument(s) ', paste(unrecog, collapse=', '))

            .Object@wsName = wsName

            # if invoked (indirectly) via a server openWs or useWs
            # method, the server will be passed in and used. if
            # invoked directly, need to create a new server instance.
            if (!is.null(serverWrap)) {
              # recycle existing server instance.
              .Object@server = serverWrap$server
            }
            else {
              # create new server instance.
              .Object@server = new('nwsServer', serverHost=serverHost, port=port)
              # now give the server a chance to do its thing.
              spaceWrap = new.env()
              spaceWrap$space = .Object
              handler = function(e) { close(.Object@server@nwsSocket); stop(e) }
              if (useUse) {
                # don't claim this space.
                tryCatch(nwsUseWs(.Object@server, wsName, spaceWrap, ...), error=handler)
              }
              else {
                # attempt to claim ownership
                tryCatch(nwsOpenWs(.Object@server, wsName, spaceWrap, ...), error=handler)
              }
            }
            .Object@cookieProtocol <- .Object@server@cookieProtocol

            .Object
          })

showNetWorkSpace <- function(object) {
    nws <- object
    server <- nws@server

    cat('\n')
    cat('NWS Host:\t', server@serverHost, ':', server@port, '\n', sep='')
    cat('Workspace Name:\t', nws@wsName, '\n', sep='')
    cat('\n')
}

setMethod('show', 'netWorkSpace', showNetWorkSpace)

setGeneric('nwsClose', function(.Object) standardGeneric('nwsClose'))
setGeneric('nwsDeclare', function(.Object, xName, mode) standardGeneric('nwsDeclare'))
setGeneric('nwsDeleteVar', function(.Object, xName) standardGeneric('nwsDeleteVar'))
setGeneric('nwsFetch', function(.Object, xName) standardGeneric('nwsFetch'))
setGeneric('nwsFetchTry', function(.Object, xName, defaultVal=NULL) standardGeneric('nwsFetchTry'))
setGeneric('nwsFind', function(.Object, xName) standardGeneric('nwsFind'))
setGeneric('nwsFindTry', function(.Object, xName, defaultVal=NULL) standardGeneric('nwsFindTry'))
setGeneric('nwsFetchFile', function(.Object, xName, fObj) standardGeneric('nwsFetchFile'))
setGeneric('nwsFetchTryFile', function(.Object, xName, fObj) standardGeneric('nwsFetchTryFile'))
setGeneric('nwsFindFile', function(.Object, xName, fObj) standardGeneric('nwsFindFile'))
setGeneric('nwsFindTryFile', function(.Object, xName, fObj) standardGeneric('nwsFindTryFile'))
setGeneric('nwsIFetch', function(.Object, xName) standardGeneric('nwsIFetch'))
setGeneric('nwsIFetchTry', function(.Object, xName, defaultVal=NULL) standardGeneric('nwsIFetchTry'))
setGeneric('nwsIFind', function(.Object, xName) standardGeneric('nwsIFind'))
setGeneric('nwsIFindTry', function(.Object, xName, defaultVal=NULL) standardGeneric('nwsIFindTry'))
setGeneric('nwsListVars', function(.Object, wsName='', showDataFrame=FALSE) standardGeneric('nwsListVars'))
setGeneric('nwsStore', function(.Object, xName, xVal) standardGeneric('nwsStore'))
setGeneric('nwsStoreFile', function(.Object, xName, fObj, n=0) standardGeneric('nwsStoreFile'))
setGeneric('nwsWsName', function(.Object) standardGeneric('nwsWsName'))
setGeneric('nwsVariable', function(.Object, xName, mode=c('fifo','lifo','multi','single'),
           env=parent.frame(), force=FALSE, quietly=FALSE) standardGeneric('nwsVariable'))
setGeneric('nwsServerObject', function(.Object) standardGeneric('nwsServerObject'))

setMethod('nwsClose', 'netWorkSpace',
          function(.Object) {
            # XXX this seems wrong
            close(.Object@server)
          })

# setGeneric('close', function(con, ...) standardGeneric('close'))
# setMethod('close','netWorkSpace', function(con) nwsClose(con))

# helper function for nwsDeclare method.
nwsDeclareInternal <- function(s, ws, xName, mode) {
  op = 'declare var'

  writeBin(charToRaw(sprintf('0004%020d%s%020d%s%020d%s%020d%s',
                             nchar(op), op,
                             nchar(ws), ws,
                             nchar(xName), xName,
                             nchar(mode), mode)), s)

  as.integer(nwsRecvN(s, 4))
}

setMethod('nwsDeclare', 'netWorkSpace',
          function(.Object, xName, mode) {
            status = nwsDeclareInternal(.Object@server@nwsSocket, .Object@wsName, xName, mode)
            if (status != 0) {
              stop('variable declaration failed')
            }
          })

setMethod('nwsDeleteVar', 'netWorkSpace',
          function(.Object, xName) {
            op = 'delete var'
            s = .Object@server@nwsSocket
            ws = .Object@wsName

            writeBin(charToRaw(sprintf('0003%020d%s%020d%s%020d%s',
                                    nchar(op), op,
                                    nchar(ws), ws,
                                    nchar(xName), xName)), s)

            status = as.integer(nwsRecvN(s, 4))
            if (status != 0) {
              stop('deleteVar failed')
            }
          })

# helper function for fetch/find methods.
nwsRetrieve <- function(cprot, s, ws, xName, op, defaultVal=NULL) {
  sn = nchar(c(op, ws, xName))
  writeBin(charToRaw(sprintf('0003%020d%s%020d%s%020d%s',
                          sn[1], op,
                          sn[2], ws,
                          sn[3], xName)), s)

  status = as.integer(nwsRecvN(s, 4))

  desc = as.integer(nwsRecvN(s, 20))
  envId = desc %/% 16777216 #(2^24)
  # if bit zero is set, then the object is not serialized
  notSerialized = desc %% 2

  if (cprot) cookie <- nwsRecvN(s, 40)

  n = as.integer(nwsRecvN(s, 20))
  sVal = nwsRecvN(s, n, rawflag=TRUE)

  if (status != 0) stop('retrieval failed')

  if (notSerialized) {
    # if bit zero and one of desc are set, it's binary data
    if (desc %% 4 - 1)
      # Return a raw vector
      sVal
    else
      # Return a character string
      rawToChar(sVal)
  }
  else if (length(sVal) > 0) {
    # Return an object
    unserialize(sVal)
  }
  else {
    # Return the defalt value
    defaultVal
  }
}

setMethod('nwsFetch', 'netWorkSpace',
          function(.Object, xName) {
            nwsRetrieve(.Object@cookieProtocol, .Object@server@nwsSocket,
                        .Object@wsName, xName, 'fetch')
          })

setMethod('nwsFetchTry', 'netWorkSpace',
          function(.Object, xName, defaultVal=NULL) {
            tryCatch(nwsRetrieve(.Object@cookieProtocol,
                                 .Object@server@nwsSocket, .Object@wsName,
                                 xName, 'fetchTry', defaultVal),
                    error=function(e) defaultVal)
          })

setMethod('nwsFind', 'netWorkSpace',
          function(.Object, xName) {
            nwsRetrieve(.Object@cookieProtocol, .Object@server@nwsSocket,
                        .Object@wsName, xName, 'find')
          })

setMethod('nwsFindTry', 'netWorkSpace',
          function(.Object, xName, defaultVal=NULL) {
            tryCatch(nwsRetrieve(.Object@cookieProtocol,
                                 .Object@server@nwsSocket, .Object@wsName,
                                 xName, 'findTry', defaultVal),
                    error=function(e) defaultVal)
          })

# helper function for fetchFile/findFile methods.
nwsRetrieveFile <- function(cprot, s, ws, xName, op, fObj) {
  if (missing(fObj)) {
    stop('no value specified for fObj argument')
  }

  if (is.character(fObj)) {
    f <- file(fObj, 'wb')
    on.exit(close(f))
  } else {
    if (!is(fObj, "file") || !isOpen(fObj, "w") || summary(fObj)$text != "binary")
      stop('fobj must be a binary mode file object opened for writing')
    f <- fObj
  }

  sn <- nchar(c(op, ws, xName))
  writeBin(charToRaw(sprintf('0003%020d%s%020d%s%020d%s',
                          sn[1], op,
                          sn[2], ws,
                          sn[3], xName)), s)

  status <- as.integer(nwsRecvN(s, 4))

  # even if failure status, read the rest of the bytes
  desc <- as.integer(nwsRecvN(s, 20))
  if (cprot) cookie <- nwsRecvN(s, 40)
  n <- as.integer(nwsRecvN(s, 20))

  blen <- 16 * 1024
  while (n > 0) {
    d <- nwsRecvN(s, min(n, blen), rawflag=TRUE)
    if (length(d) == 0) stop('NWS server connection dropped')
    writeBin(d, f)
    n <- n - length(d)
  }

  if (status != 0) stop('retrieval failed')
  TRUE
}

setMethod('nwsFetchFile', 'netWorkSpace',
          function(.Object, xName, fObj) {
            nwsRetrieveFile(.Object@cookieProtocol, .Object@server@nwsSocket,
                            .Object@wsName, xName, 'fetch', fObj)
          })

setMethod('nwsFetchTryFile', 'netWorkSpace',
          function(.Object, xName, fObj) {
            tryCatch({
                nwsRetrieveFile(.Object@cookieProtocol, .Object@server@nwsSocket,
                                .Object@wsName, xName, 'fetchTry', fObj)
              }, error=function(e) {
                if (e$message == 'retrieval failed') {
                  FALSE
                }
                else {
                  stop(e$message)
                }
              })
          })

setMethod('nwsFindFile', 'netWorkSpace',
          function(.Object, xName, fObj) {
            nwsRetrieveFile(.Object@cookieProtocol, .Object@server@nwsSocket,
                            .Object@wsName, xName, 'find', fObj)
          })

setMethod('nwsFindTryFile', 'netWorkSpace',
          function(.Object, xName, fObj) {
            tryCatch({
                nwsRetrieveFile(.Object@cookieProtocol, .Object@server@nwsSocket,
                                .Object@wsName, xName, 'findTry', fObj)
              }, error=function(e) {
                if (e$message == 'retrieval failed') {
                  FALSE
                }
                else {
                  stop(e$message)
                }
              })
          })

setMethod('nwsIFetch', 'netWorkSpace',
          function(.Object, xName) {
            nwsValueIterator(.Object, xName, 'ifetch', NULL)
          })

setMethod('nwsIFetchTry', 'netWorkSpace',
          function(.Object, xName, defaultVal=NULL) {
            nwsValueIterator(.Object, xName, 'ifetchTry', defaultVal)
          })

setMethod('nwsIFind', 'netWorkSpace',
          function(.Object, xName) {
            nwsValueIterator(.Object, xName, 'ifind', NULL)
          })

setMethod('nwsIFindTry', 'netWorkSpace',
          function(.Object, xName, defaultVal=NULL) {
            nwsValueIterator(.Object, xName, 'ifindTry', defaultVal)
          })

# to see list output clearly use: write(nwsList...(), stdout())
setMethod('nwsListVars', 'netWorkSpace',
          function(.Object, wsName='', showDataFrame=FALSE) {
            op = 'list vars'
            s = .Object@server@nwsSocket
            if (wsName == '') wsName = .Object@wsName

            writeBin(charToRaw(sprintf('0002%020d%s%020d%s',
                                    nchar(op), op,
                                    nchar(wsName), wsName)), s)

            # status, unused at the moment
            status = as.integer(nwsRecvN(s, 4))
            desc = nwsRecvN(s, 20)
            if (.Object@cookieProtocol)
              cookie <- nwsRecvN(s, 40)

            ret <- nwsRecvN(s, as.integer(nwsRecvN(s, 20)))
            if (showDataFrame==FALSE)
              ret
            else {
              ## convert response into an R data frame
              ret <- unlist(strsplit(ret, "\n"))
              retval <- list()
              fields <- list()

              i = 1
              while (i<=length(ret)) {
                line <- unlist(strsplit(ret[i], "\t"))

                # convert each field to correct type
                fields[1] = line[1]
                fields[2] = as.integer(line[2])
                fields[3] = as.integer(line[3])
                fields[4] = as.integer(line[4])
                fields[5] = line[5]
                retval = c(retval, list(fields))
                i = i+1
              }

              if (length(retval)>0) {
                names(retval) <- seq(along=retval)
                retval <- do.call(rbind, retval)
                colnames(retval) <-
                c("Variable", "NumValues", "NumFetchers", "NumFinders", "Mode")
              }

              retval <- data.frame(retval)
              retval
            }

          })

# helper function for store method
nwsStoreInternal <- function(s, ws, xName, xVal) {
  op = 'store'

  desc = nwsRFP # R Fingerprint

  if (missing(xVal)) {
    stop('no value specified for xVal argument')
  }

  if (is.raw(xVal)) {
    desc = desc + 3
  }
  else if (!is.character(xVal) || (length(xVal) != 1)) {
    xVal = serialize(xVal, ascii=FALSE, connection=NULL)
    # serialize returns a raw vector as of R 2.4 
    if (is.character(xVal)) xVal = charToRaw(xVal)
  }
  else {
    xVal = charToRaw(xVal)
    desc = desc + 1 # in other systems, we use a manifest constant and a bit or here... .
  }
  descTxt = sprintf('%020i', desc) # would prefer to use unsigned here.

  sn = nchar(c(op, ws, xName, descTxt))
  xLen = length(xVal)
  writeBin(c(charToRaw(sprintf('0005%020d%s%020d%s%020d%s%020d%s%020d',
                             sn[1], op,
                             sn[2], ws,
                             sn[3], xName,
                             sn[4], descTxt,
                             xLen)), xVal), s)

  # status, barely used at the moment.
  status = as.integer(nwsRecvN(s, 4))

  if (status != 0) {
    stop('store failed')
  }
}

setMethod('nwsStore', 'netWorkSpace',
          function(.Object, xName, xVal) {
            nwsStoreInternal(.Object@server@nwsSocket, .Object@wsName, xName, xVal)
          })

setMethod('nwsStoreFile', 'netWorkSpace',
          function(.Object, xName, fObj, n=0) {
            ws <- .Object@wsName
            s <- .Object@server@nwsSocket
            op <- 'store'

            desc <- nwsRFP + 3  # R Fingerprint and raw data

            if (missing(fObj)) {
              stop('no value specified for fObj argument')
            }

            # if fObj is a character string, handle it specially
            if (is.character(fObj)) {
              f <- file(fObj, 'rb')
              on.exit(close(f))
            } else {
              if (!is(fObj, "file") || !isOpen(fObj, "r") || summary(fObj)$text != "binary")
                stop('fobj must be a binary mode file object opened for reading')
              f <- fObj
            }
      
            fsize <- file.info(summary(f)$description)$size
            fpos <- seek(f)
            fbytes <- fsize - fpos
            n <- if (n <= 0) fbytes else min(n, fbytes)
            if (n <= 0) return(FALSE)

            descTxt <- sprintf('%020i', desc) # would prefer to use unsigned here.

            sn <- nchar(c(op, ws, xName, descTxt))
            writeBin(charToRaw(sprintf('0005%020d%s%020d%s%020d%s%020d%s%020d',
                                       sn[1], op,
                                       sn[2], ws,
                                       sn[3], xName,
                                       sn[4], descTxt,
                                       n)), s)

            blen <- 16 * 1024
            while (n > 0) {
              d <- readBin(f, what='raw', n=min(blen, n))
              dlen <- length(d)
              if (dlen <= 0)
                break
              writeBin(d, s)
              n <- n - dlen
            }

            if (n > 0) {
              # I don't thing this should ever happen unless the file
              # size computation is incorrect, but I really don't want
              # to corrupt the server connection
              warning('unable to read all the data in file ',
                      summary(f)$description, ' [size: ', fsize,
                      'bytes]: padding value in workspace variable')
              blen <- 1024
              buffer <- raw(blen)
              while (n > 0) {
                if (blen <= n) {
                  writeBin(buffer, s)
                  n <- n - dlen
                } else {
                  writeBin(raw(n), s)
                  break
                }
              }
            }

            # status, barely used at the moment.
            status <- as.integer(nwsRecvN(s, 4))

            if (status != 0) {
              stop('store file failed')
            }
            TRUE
          })

setMethod('nwsWsName', 'netWorkSpace', function(.Object) {.Object@wsName})

setMethod('nwsVariable', 'netWorkSpace',
          function(.Object, xName, mode=c('fifo','lifo','multi','single'),
                   env=parent.frame(), force=FALSE, quietly=FALSE) {
            missingMode = missing(mode)
            mode <- match.arg(mode)

            # be careful, because 'exists' will cause an active binding function
            # to be called, which is a side effect that we don't want
            if (force ||
                (!tryCatch(bindingIsActive(xName, env), error=function(...) FALSE) &&
                 !exists(xName, envir=env, inherits=FALSE))) {
              s = .Object@server@nwsSocket
              ws = .Object@wsName
              if (missingMode) {
                mlist = c(mode, 'fifo', 'lifo', 'multi', 'single')
                mode = NA
                for (m in mlist) {
                  if (nwsDeclareInternal(s, ws, xName, m) == 0) {
                    mode = m
                    break
                  }
                }
                if (is.na(mode))
                  stop('unable to declare variable')
              } else {
                if (nwsDeclareInternal(s, ws, xName, mode) != 0)
                  stop('variable declaration failed')
              }

              if (identical(mode, 'single')) {
                mf <- function(val)
                  if (missing(val))
                    nwsRetrieve(.Object@cookieProtocol, s, ws, xName, 'find')
                  else
                    nwsStoreInternal(s, ws, xName, val)
              } else {
                mf <- function(val)
                  if (missing(val))
                    nwsRetrieve(.Object@cookieProtocol, s, ws, xName, 'fetch')
                  else
                    nwsStoreInternal(s, ws, xName, val)
              }

              t <- makeActiveBinding(xName, mf, env)
            } else {
              if (! quietly)
                warning('not overwriting previous binding for ', xName)
            }
          })

setMethod('nwsServerObject', 'netWorkSpace', function(.Object) .Object@server)

# helper function for ifetch/ifind methods.
nwsIRetrieve <- function(s, ws, xName, op, varId, valIndex) {
  sn = nchar(c(op, ws, xName))
  writeBin(charToRaw(sprintf('0005%020d%s%020d%s%020d%s%020d%-20.20s%020d%020d',
                          sn[1], op,
                          sn[2], ws,
                          sn[3], xName,
                          20, varId,
                          20, valIndex)), s)

  status = as.integer(nwsRecvN(s, 4))

  desc = as.integer(nwsRecvN(s, 20))
  envId = desc %/% 16777216 #(2^24)
  # if bit zero is set, then the object is not serialized
  notSerialized = desc %% 2

  # cookie protocol is assumed at this point
  varId = nwsRecvN(s, 20)
  valIndex = as.integer(nwsRecvN(s, 20))

  n = as.integer(nwsRecvN(s, 20))
  sVal = nwsRecvN(s, n, rawflag=TRUE)

  if (notSerialized) {
    # if bit zero and one of desc are set, it's binary data
    if (desc %% 4 - 1)
      # Return a raw vector
      list(status=status, sVal=sVal, varId=varId, valIndex=valIndex)
    else
      # Return a character string
      list(status=status, sVal=rawToChar(sVal), varId=varId, valIndex=valIndex)
  }
  else if (length(sVal) > 0) {
    list(status=status, sVal=unserialize(sVal), varId=varId, valIndex=valIndex)
  }
  else {
    stop('StopIteration')
  }
}

# helper function to return a closure that acts as an iterator
nwsValueIterator <- function(.Object, xName, op, defaultVal) {
  if (!.Object@cookieProtocol)
    stop('NWS server does not support iterated operations')
  if (!is.character(xName))
    stop('variable name must be a string')

  # initial state of the closure
  varId <- ''
  valIndex <- 0

  function() {
    defval <- list(status=0, varId=varId, valIndex=valIndex, sVal=defaultVal)

    r <- tryCatch({
          nwsIRetrieve(.Object@server@nwsSocket, .Object@wsName,
                       xName, op, varId, valIndex)
        }, error=function(e) defval)

    varId <<- r$varId
    valIndex <<- r$valIndex

    if (r$status != 0) stop('retrieval failed')
    r$sVal
  }
}
