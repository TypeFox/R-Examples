#############################################################################
#
# Copyright 2015 Axibase Corporation or its affiliates. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License").
# You may not use this file except in compliance with the License.
# A copy of the License is located at
#
# https://www.axibase.com/atsd/axibase-apache-2.0.pdf
#
# or in the "license" file accompanying this file. This file is distributed
# on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
# express or implied. See the License for the specific language governing
# permissions and limitations under the License.
#
#############################################################################
#
#' Show connection parameters.
#'
#' @description
#' The function shows the current values of the connection parameters
#' \code{url, user, password, verify} and \code{encryption}.
#' They are used to arrange a connection with ATSD.
#' @seealso
#' You could change the connection parameters with the \code{\link{set_connection}}
#' function and save that changes to the configuration file with
#' the \code{\link{save_connection}} function.
#' @export
show_connection <- function() {
  cat(" url = ", get("url", envir = atsdEnv))
  cat("\n user = ", get("user", envir = atsdEnv))
  cat("\n password = ", get("password", envir = atsdEnv))
  verify <- get("verify", envir = atsdEnv)
  if (is.na(verify)) {
    cat("\n verify = NA")
  } else if (verify) {
    cat("\n verify = yes")
  } else {
    cat("\n verify = no")
  }
  cat("\n encryption = ", get("encryption", envir = atsdEnv))
}
