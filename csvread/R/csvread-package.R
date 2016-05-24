#-------------------------------------------------------------------------------
#
# Package csvread 
#
# General description 
# 
# Sergei Izrailev, 2011-2015
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# Copyright 2015 Jabiru Ventures LLC
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
#' Package \code{csvread} contains a fast specialized CSV and other delimited 
#' file loader, and a basic 64-bit integer class to aid in reading 64-bit 
#' integer values.
#' 
#' @name csvread
#' @aliases csvread
#' @docType package
#' @title Fast Specialized CSV File Loader. 
#' @author Sergei Izrailev
#' @section Maintainer: Sergei Izrailev 
#' @section Copyright: Copyright (C) Collective, Inc.; with portions Copyright (C) Jabiru Ventures LLC
#' @section License: Apache License, Version 2.0, 
#'    available at http://www.apache.org/licenses/LICENSE-2.0
#' @section URL: http://github.com/jabiru/csvread
#' @section Installation from github: 
#' \code{devtools::install_github("jabiru/csvread")}
#' @keywords csv delimited file read.csv bigint 64-bit integer64
#' @useDynLib csvread
#' @exportPattern "*"
#' @import methods
#' 
# The next and last line should be the word 'NULL'.
NULL
