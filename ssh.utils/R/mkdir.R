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

#' Creates a remote directory with the specified group ownership and permissions.
#' 
#' If the directory already exists, attempts to set the group ownership to the 
#' \code{user.group}. The allowed group permissions are one of 
#' \code{c("g+rwx", "g+rx", "go-w", "go-rwx")}, or \code{"-"}. The value 
#' \code{"-"} means "don't change permissions". 
#' @param path Directory path. If using \code{remote}, this should be a full path or 
#'             a path relative to the user's home directory.
#' @param user.group The user group. If NULL, the default group is used.
#' @param remote Remote machine specification for ssh, in format such as \code{user@@server} that does not 
#'        require interactive password entry. For local execution, pass an empty string "" (default).
#' @param permissions The group permissions on the directory. Default is 'rwx'.
#' @rdname mkdir.remote
#' @note This may not work on Windows.
# COMPATIBILITY WARNING: WINDOWS 
mkdir.remote <- function(path, user.group = NULL, remote = "", 
      permissions = c("g+rwx", "g+rx", "go-w", "go-rwx", "-")) 
{
   permissions <- match.arg(permissions)
   if (!file.exists.remote(path, remote = remote)) 
   {
      run.remote(paste("mkdir -p", path), remote = remote)
   }
   if (!is.null(user.group))
   {
      run.remote(paste("chgrp -R", user.group, path), remote = remote)
   }
   if (permissions != "-")
   {
      run.remote(paste("chmod ", permissions, " ", path, sep = ""), remote = remote)         
   }
}

#-------------------------------------------------------------------------------
