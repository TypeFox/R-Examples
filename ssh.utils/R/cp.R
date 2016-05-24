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

#' A wrapper around the scp shell command that handles local/remote files and allows 
#' copying between remote hosts via the local machine.
#' 
#' @param remote.src Remote machine for the source file in the format \code{user@@machine} or an empty string for local.
#' @param path.src Path of the source file.
#' @param remote.dest Remote machine for the destination file in the format \code{user@@machine} or an empty string for local.
#' @param path.dest Path for the source file; can be a directory.
#' @param verbose Prints elapsed time if TRUE
#' @param via.local Copies the file via the local machine. Useful when two remote machines can't talk to each other directly.
#' @param local.temp.dir When copying via local machine, the directory to use as scratch space. 
#' @name cp.remote
#' @aliases cp.remote
#' @rdname cp.remote
#' @title scp wrapper
#' @examples
#' \dontrun{
#' ## Copy file myfile.csv from the home directory on the remote server to
#' ## the local working directory.
#' 
#' ## on remote server in bash shell:
#' # cat myfile.csv
#' # [me@@myserver ~]$ cat myfile.csv
#' # "val","ts"
#' # 1,
#' # 2,
#' # 3,
#' # 4,
#' # 5,
#' # 6,
#' # 7,
#' # 8,
#' # 9,
#' # 10,
#' 
#' ## on local server in R:
#' cp.remote(remote.src = "me@@myserver", path.src = "~/myfile.csv", 
#'           remote.dest = "", path.dest = getwd(), verbose = TRUE)
#' # [1] "Elapsed: 1.672 sec"
#' df <- read.csv("myfile.csv")
#' df
#' #    val ts
#' # 1    1 NA
#' # 2    2 NA
#' # 3    3 NA
#' # 4    4 NA
#' # 5    5 NA
#' # 6    6 NA
#' # 7    7 NA
#' # 8    8 NA
#' # 9    9 NA
#' # 10  10 NA
#' }
cp.remote <- function(remote.src, path.src, remote.dest, path.dest, verbose = FALSE, 
      via.local = FALSE, local.temp.dir = tempdir())
{   
   if (remote.src == "" || is.na(remote.src)) remote.src = NULL;
   path.src <- paste(c(remote.src, path.src), collapse = ":")
   
   if (remote.dest == "" || is.na(remote.dest)) remote.dest = NULL;
   path.dest <- paste(c(remote.dest, path.dest), collapse = ":")
   
   if (via.local)
   {
      path.tmp <- tempfile(pattern = "tmp_scp_", tmpdir = local.temp.dir) 
      res <- run.remote(paste("scp", path.src, path.tmp), remote="")
      if (res$cmd.error)
      {
         stop(paste("scp failed:", res$warn.msg))
      }
      if (verbose) print(paste("Elapsed:", attr(res$cmd.out, "elapsed.time"), "sec"))
      path.src <- path.tmp
   }
   res <- run.remote(paste("scp", path.src, path.dest), remote="")
   if (res$cmd.error)
   {
      stop(paste("scp failed:", res$warn.msg))
   }
   if (verbose) print(paste("Elapsed:", attr(res$cmd.out, "elapsed.time"), "sec"))   
   
   if (via.local) run.remote(paste("rm -f", path.tmp), remote="")
}

#-------------------------------------------------------------------------------
