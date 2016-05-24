#-------------------------------------------------------------------------------
#
# Package ssh.utils 
#
# Implementation. 
# 
# Sergei Izrailev, 2011-2014
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------

#' Checks for processes running on a local or remote machine.
#' 
#' One of the use cases for this function is to ensure that an R process is 
#' already running and not start another one accidentally.
#' @param grep.string  String(s) to check for in \code{ps}. If a vector, runs a chain of piped grep commands for each string.
#' @param remote Remote machine specification for ssh, in format such as \code{user@@server} that does not 
#'        require interactive password entry. For local execution, pass an empty string "" (default).
#' @param stop.if.any  Stop if any of \code{grep.string} is running
#' @param stop.if.none Stop if none of \code{grep.string} is running
#' @param count.self When \code{FALSE}, excludes the calling process name from the count, if it gets matched. 
#' @param ps.options   Gives the ability to run different options to ps.
#' @seealso \code{run.remote}
#' @rdname ps.grep.remote
#' @note This may not work on Windows.
#' @examples 
#' \dontrun{
#' # Check if Eclipse is running.
#' ps.grep.remote("Eclipse", remote = "")
#' # [1] TRUE
#' }
ps.grep.remote <- function(
      grep.string             # string(s) to check for in 'ps'
      , remote                # remote location to look for script
      , stop.if.any = FALSE   # stop if any of 'grep.string' running
      , stop.if.none = FALSE  # stop if none of 'grep.string' running
      , count.self = FALSE    # if TRUE, will count all processes including itself; could be useful if remote is different
      , ps.options = "aux"
) 
{
   if (length(grep.string) == 0) stop("ps.grep.remote: grep.string is empty")
   
   # Assemble command
   greps <- unlist(lapply(grep.string, function(s) { paste("egrep '", s, "'", sep = "") } ))
   pipes <- paste(greps, collapse = " | ")
   cmd <- paste("ps ", ps.options, " | ", pipes, " | grep -v 'egrep'", sep = "")
   
   # Execute ps remotely   
   res <- run.remote(cmd, remote = remote)$cmd.out
   
   # Filter again for grep.string
   res <- res[str_detect(res, grep.string)]
   
   # Check number of rows
   running <- length(res) > 0
   
   if (running && !count.self)
   {
      # check if what's running is the process that called this function (self)
      script.call <- paste(commandArgs(), collapse=" ")
      res.match <- sapply(res, function(s) length(grep(script.call, s, ignore.case = T, perl = T, value = T)) > 0)
      # either we didn't match (e.g. on remote host) or we matched more than one
      running <- !(sum(res.match) == 1)
   }
   
   # Stop if the 'grep.string' is running   
   if (stop.if.any && running)
   {
      stop(paste("Found ", grep.string, " in ps aux on remote server '", remote, "'; stopping\n\n",
                  paste(res, collapse = "\n"),
                  sep = ""))
   }
   
   # Stop if the 'grep.string' is not running
   if (stop.if.none && !running)
   {
         stop(paste("Did not find ", grep.string, " in ps aux on remote server '", remote, "'; stopping", sep = ""))      
   }
   
   return(running)   
}

#------------------------------------------------------------------------------
