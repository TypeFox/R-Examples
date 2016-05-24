#
#   This work was created by the National Center for Ecological Analysis and Synthesis.
#
#     Copyright 2015 Regents of the University of California
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
#' Determine whether an externalptr object is NULL.
#' 
#' The pointer is treated as an externalptr and checked via a call in C to see if it is NULL.
#' @param pointer externalptr to be checked for NULL value
#' @return logical TRUE if the pointer is NULL, otherwise FALSE
#' @export
is.null.externalptr <- function(pointer) {
    if (class(pointer)!="externalptr") stop("pointer is not an externalptr.")
    rv <- .Call("isnull", pointer)
    return(rv)
}
