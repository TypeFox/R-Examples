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

#' Checks if a local or remote file exists.
#' 
#' A wrapper around a bash script. Works with local files too if \code{remote=""}.
#' @param file File path.
#' @param remote Remote machine specification for ssh, in format such as \code{user@@server} that does not 
#'        require interactive password entry. For local execution, pass an empty string "" (default).
#' @return \code{TRUE} or \code{FALSE} indicating whether the file exists.
#' @rdname file.exists.remote
#' @examples 
#' \dontrun{
#' file.exists.remote("~/myfile.csv", remote = "me@@myserver")
#' # [1] TRUE
#' }
file.exists.remote <- function(file, remote = "")
{
   cmd <- paste("if [ -e ", file, " ] ; then echo TRUE; else echo FALSE; fi ", sep="")
   res <- run.remote(cmd, remote)
   if (res$cmd.error)
   {
      stop(paste("file.exists.remote: ERROR", res$cmd.out, res$warn.msg, sep="\n"))
   }
   as.logical(res$cmd.out)
}

#-------------------------------------------------------------------------------
