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

#' Measure the resident memory usage of a process.
#' 
#' Returns the memory usage in KB of a process with the specified process id. By 
#' default, returns the memory usage of the current R process. This can be used
#' to measure and log the memory usage of the R process during script execution.
#' @param pid Process ID (default is the current process id).
#' @return The resident memory usage in KB. 
#' @examples
#' \dontrun{
#' mem.usage()
#' # [1] 37268
#' }
#' @rdname mem.usage
mem.usage <- function(pid = Sys.getpid()) 
{ 
   df <- read.delim(pipe("ps axo pid,rss"), sep = ""); 
   df[df$PID == pid, "RSS"] 
}

#-------------------------------------------------------------------------------
