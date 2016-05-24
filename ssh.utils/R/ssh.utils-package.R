#-------------------------------------------------------------------------------
#
# Package ssh.utils 
#
# General description 
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
 
#' Package \code{ssh.utils} provides utility functions for system command
#' execution, both locally and remotely using ssh/scp. The command 
#' output is captured and provided to the caller. This functionality is
#' intended to streamline calling shell commands from R, retrieving and 
#' using their output, while instrumenting the calls with appropriate
#' error handling. NOTE: this first version is limited to unix with local 
#' and remote systems running bash as the default shell. 
#' 
#' @name ssh.utils
#' @aliases ssh.utils
#' @seealso \code{\link{run.remote}}, \code{\link{cp.remote}}, 
#'          \code{\link{file.exists.remote}}, \code{\link{mkdir.remote}},
#'          \code{\link{ps.grep.remote}}, \code{\link{mem.usage}} 
#' @docType package
#' @title Local and remote system commands with output and error capture. 
#' @author Sergei Izrailev
#' @section OS_type: unix 
#' @section Maintainer: Sergei Izrailev 
#' @section Copyright: Copyright (C) Collective, Inc.
#' @section License: Apache License, Version 2.0, 
#'    available at http://www.apache.org/licenses/LICENSE-2.0
#' @section URL: http://github.com/collectivemedia/ssh.utils
#' @section Installation from github: 
#' \code{devtools::install_github("collectivemedia/ssh.utils")}
#' @keywords ssh scp remote shell bash "capture output" "system command"
#' @exportPattern "*"
#' @import stringr 
#' 
# The next and last line should be the word 'NULL'.
NULL
