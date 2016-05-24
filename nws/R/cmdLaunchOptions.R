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

sshcmd <- function(host, options) {
  basicArgs <- if (!is.null(options$user))
                 c('ssh', '-f', '-x', '-l', options$user, host)
               else
                 c('ssh', '-f', '-x', host)

  wrapper <- file.path(options$wrapperDir, options$workerWrapper)
  if (file.access(wrapper) == 0) {
    if (identical(grep('\\.py$', wrapper, ignore.case=TRUE), as.integer(1))) {
      if (!is.null(options$python))
        c(options$python, wrapper, basicArgs)
      else
        c('python', wrapper, basicArgs)
    }
    else {
      c(wrapper, basicArgs)
    }
  }
  else {
    basicArgs
  }
}

sshforwardcmd <- function(host, options) {
  if (is.null(options$nwsHostRemote))
    stop('must use the nwsHostRemote option with sshforwardcmd')

  r <- if (nchar(options$nwsHostRemote) > 0)
         sprintf('%s:%d:%s:%d', options$nwsHostRemote, options$nwsPortRemote,
                 options$nwsHost, options$nwsPort)
       else
         sprintf('%d:%s:%d', options$nwsPortRemote,
                 options$nwsHost, options$nwsPort)

  basicArgs <- if (!is.null(options$user))
                 c('ssh', '-f', '-x', '-R', r, '-l', options$user, host)
               else
                 c('ssh', '-f', '-x', '-R', r, host)

  wrapper <- file.path(options$wrapperDir, options$workerWrapper)
  if (file.access(wrapper) == 0) {
    if (identical(grep('\\.py$', wrapper, ignore.case=TRUE), as.integer(1))) {
      if (!is.null(options$python))
        c(options$python, wrapper, basicArgs)
      else
        c('python', wrapper, basicArgs)
    }
    else {
      c(wrapper, basicArgs)
    }
  }
  else {
    basicArgs
  }
}

rshcmd <- function(host, options) {
  basicArgs <- if (!is.null(options$user))
                 c('rsh', host, '-l', options$user, '-n')
               else
                 c('rsh', host, '-n')

  wrapper <- file.path(options$wrapperDir, 'BackgroundLaunch.py')
  if (!is.null(options$python))
    c(options$python, wrapper, basicArgs)
  else
    c('python', wrapper, basicArgs)
}

lsfcmd <- function(host, options) {
  'bsub'
}

ccscmd <- function(host, options) {
  c("job", "submit", "/exclusive:false")
}

rwincmd <- function(host, options) {
  # Note: Execution of cscript (done locally) must use simple quoting.
  # However, the remote command (the one executed via rwin.vbs) must
  # be done with MSC quoting, since the remote command is presumed to
  # be the Python interpretter.  Therefore, we set options$simpleQuote
  # to true, but we don't use the rwin.vbs "-s" option.
  options$simpleQuote <- TRUE

  wrapper <- file.path(options$wrapperDir, 'rwin.vbs')
  if (is.null(options$passwd)) {
    c('cscript', '//nologo', wrapper, host, '--')
  }
  else {
    user <- if (is.null(options$user)) Sys.info()[['login']] else options$user
    c('cscript', '//nologo', wrapper, host, '-l', user, '-p', options$passwd, '--')
  }
}

envcmd <- function(host, envVars, options) {
  c('env', envVars, file.path(options$scriptDir, options$scriptName))
}

scriptcmd <- function(host, envVars, options) {
  if (!is.null(options$python))
    c(options$python, file.path(options$scriptDir, options$scriptName), envVars)
  else
    c('python', file.path(options$scriptDir, options$scriptName), envVars)
}
