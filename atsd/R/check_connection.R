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
#' @keywords internal
check_connection <- function() {
  url <- get("url", envir = atsdEnv)
  user <- get("user", envir = atsdEnv)
  password <- get("password", envir = atsdEnv)
  cert <- get("verify", envir = atsdEnv)
  encrypt <- get("encryption", envir = atsdEnv)
  if (is.na(url) || is.na(user) || is.na(password)) {
    cat("\n Wrong connection, current values of connection parameters are:\n")
    show_connection()
    cat("\n Use set_connection() and save_connection() functions to set up connection.")
    return(FALSE)
  } else if (substr(url, 1, 5) == "https" &&
               (!(cert %in% c(TRUE, FALSE)) ||
                !(encrypt %in% c("default", "tls1", "ssl2", "ssl3")))) {
    cat("\n Wrong connection, current values of connection parameters are:\n")
    show_connection()
    cat("\n Please provide valid 'verify' and 'encription' parameters to use https protocol.")
    cat('\n verify = "yes" ensures validation of ATSD SSL certificate and verify = "no"  suppresses the validation.')
    cat('\n The encription specifies cryptographic protocol used by ATSD https server.')
    cat('\n Possible values of encription are: "default", "ssl2", "ssl3", and "tls1"') 
    cat("\n Use save_connection() and set_connection() functions to set up connection parameters.")
    return(FALSE)
  } else {
    return(TRUE)
  }
}
