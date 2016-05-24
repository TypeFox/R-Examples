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
#' Write connection parameters to configuration file.
#'
#' @description
#'     The function writes the connection parameters into configuration file.
#'
#' @inheritParams set_connection
#' 
#' @details
#'     If you call \code{save_connection()} without arguments, 
#'     then the current values of the connection parameters  (including NAs)
#'     will be written to the configuration file.
#'     If you provide some arguments, they will be written 
#'     into the configuration file. If configuration file is absent 
#'     it will be created in the atsd package folder.
#'  
#' @seealso
#'     For more information about the configuration file
#'     view the package vignette:
#'     \code{browseVignettes(package = "atsd")}.
#'           
#' @examples
#' # Write the current values of the connection parameters to the configuration file
#' save_connection()
#' 
#' # Write the user name and the password to the configuration file
#' save_connection(user = "user001", password = "123456")
#' 
#' # Write all parameters nedeed for https connection to the configuration file
#' save_connection(url = "https://my.company.com:8443", user = "user001", password = "123456", 
#'                 verify = "no", encryption = "ssl3")
#' 
#' @export
save_connection <- function(url = NA, user = NA, password = NA, 
                            verify = NA, encryption = NA) {
  
  if (all(is.na(c(url, user, password, verify, encryption)))) {
    url <- get("url", envir = atsdEnv)
    user <- get("user", envir = atsdEnv)
    password <- get("password", envir = atsdEnv)
    verify <- get("verify", envir = atsdEnv)
    if (!is.na(verify)) {
      if (verify) {
        verify <- "yes"
      } else {
        verify <- "no"
      }
    }
    encryption <- get("encryption", envir = atsdEnv)
  } else {
    connection_data <- get_connection_data()
    if (is.na(url)) {
      url <- connection_data[1]
    }
    if (is.na(user)) {
      user <- connection_data[2]
    }
    if (is.na(password)) {
      password <- connection_data[3]
    }
    if (is.na(verify)) {
      verify <- connection_data[4]
    }
    if (is.na(encryption)) {
      encryption <- connection_data[5]
    }
  }

  if (is.na(url)) {url <- ""}
  if (is.na(user)) {user <- ""}
  if (is.na(password)) {password <- ""}
  if (is.na(verify)) {verify <- ""}
  if (is.na(encryption)) {encryption <- ""}
  
  l <- character(16)
  l[1] <- "# the url of ATSD including port number"
  l[2] <- paste0("url=", url)
  l[3] <- ""
  l[4] <- "# the user name"
  l[5] <- paste0('user=', user)
  l[6] <- ""
  l[7] <- "# the user's password"
  l[8] <- paste0("password=", password)
  l[9] <- ""
  l[10] <- "# validate ATSD SSL certificate: yes, no"
  l[11] <- paste0("verify=", tolower(verify))
  l[12] <- ""
  l[13] <- "# cryptographic protocol used by ATSD https server:"
  l[14] <- "# default, ssl2, ssl3, tls1"
  l[15] <- paste0("encryption=", tolower(encryption))
  l[16] <- ""

  if (!file.exists(system.file("connection.config", package = "atsd"))) {
    if (!file.create(paste0(system.file(package = "atsd"), "/connection.config"))) {
      message("Can not create new configuration file.")
      invisible()
    } else {
      message("New configuration file is created.")
    }
  }
  fileConn <- file(system.file("connection.config", package = "atsd"))
  on.exit(close(fileConn))
  writeLines(l, fileConn)
}
