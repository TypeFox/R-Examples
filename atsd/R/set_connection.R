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
#' Set up parameters of a connection with ATSD.
#'
#' @description
#'     The function overrides the connection parameters
#'     for the duration of the current R session without changing the configuration file.
#'     
#' @param url
#'     Optional string argument. The url of ATSD with the port number.
#'     
#' @param user
#'     Optional string argument. The user name.
#'     
#' @param password
#'     Optional string argument. The user's password.
#'     
#' @param  verify
#'      Optional string argument -- "yes" or "no". 
#'      \code{verify = "yes"}  ensures validation of ATSD SSL certificate and 
#'      \code{verify = "no"}  suppresses the validation 
#'      (applicable in the case of 'https' protocol).
#'      
#' @param  encryption
#'      Optional string argument.
#'      Cryptographic protocol used by ATSD https server. 
#'      Possible values are: "default", "ssl2", "ssl3", and "tls1" 
#'      (In most cases, use "ssl3" or "tls1".)
#'      
#' @param  file
#'      Optional string argument.
#'      The absolute path to the file from which the connection parameters could be read. 
#'      The file should be formatted as the package configuration file, 
#'      see the Details section below.
#'      
#' @details
#'     The function overrides the connection parameters for the duration 
#'     of the current R session without changing the configuration file. 
#'     If called without arguments the function sets the connection parameters 
#'     from the configuration file. If the file  argument is provided the function use it. 
#'     In both cases the current values of the parameters became the same as in the file. 
#'     The file should be a plain text file formatted as the following:
#'      
#'      \code{# the url of ATSD including port number}
#'      \cr
#'      \code{url=http://host_name:port_number}
#'      \cr
#'      \code{# the user name}
#'      \cr
#'      \code{user=atsd_user_name}
#'      \cr
#'      \code{# the user's password}
#'      \cr
#'      \code{password=atsd_user_password}
#'      \cr
#'      \code{# validate ATSD SSL certificate: yes, no}
#'      \cr
#'      \code{verify=no}
#'      \cr
#'      \code{# cryptographic protocol used by ATSD https server:}
#'      \cr
#'      \code{# default, ssl2, ssl3, tls1}
#'      \cr
#'      \code{encryption=ssl3}
#'      
#'     In case the \code{file}  argument is not provided, but some of other arguments are specified, 
#'     the only specified parameters will be changed.
#' 
#' @seealso   
#'     To see the current values of the connection parameters use
#'     the \code{\link{show_connection}} function.
#'     To change the configuration file use
#'     the \code{\link{save_connection}} function.
#'     
#' @examples
#' # Modify the user
#' set_connection(user = "user001")
#' 
#' # Modify the cryptographic protocol
#' set_connection(encryption = "tls1")
#' 
#' # Set up url, user and password
#' set_connection(url = "http://my.company.com:8088", user = "user001", password = "123456")
#'     
#' # Set up parameters of https connection
#' set_connection(url = "https://my.company.com:8443", user = "user001", password = "123456", 
#'                verify = "no", encryption = "ssl3")
#' 
#' \dontrun{
#' # Set up parameters from a file
#' set_connection(file = "/home/user001/atsd_https_connection.txt")
#'  
#' # Set up parameters from the configuration file
#' set_connection()
#' }
#'     
#' @export
set_connection <- function(url = NA, user = NA, password = NA, 
                           verify = NA, encryption = NA, file = NA) {
  from <- character()
  rewrite_na <- TRUE
  if (!is.na(file)) {
    connection_data <- get_connection_data(file)
    from <- paste0("file ", file)
  } else if (!all(is.na(c(url, user, password, verify, encryption)))) {
    connection_data <- c(url, user, password, verify, encryption, "true")
    from <- "arguments"
    rewrite_na <- FALSE
  } else {
    connection_data <- get_connection_data()
    from <- paste0("configuration file ",
                   system.file("connection.config", package = "atsd"))
  }
  if (!is.na(connection_data[1]) || rewrite_na) {
    assign("url", connection_data[1], envir = atsdEnv)
  }
  if (!is.na(connection_data[2]) || rewrite_na) {
    assign("user", connection_data[2], envir = atsdEnv)
  }
  if (!is.na(connection_data[3]) || rewrite_na) {
    assign("password", connection_data[3], envir = atsdEnv)
  }
  if (!is.na(connection_data[4])) {
    if (is.character(connection_data[4]) && (tolower(connection_data[4]) %in% c("yes", "true"))) {
      assign("verify", TRUE, envir = atsdEnv)
    } else {
      assign("verify", FALSE, envir = atsdEnv)
    }
  } else if (rewrite_na) {
    assign("verify", NA, envir = atsdEnv)
  }
  if (!is.na(connection_data[5]) || rewrite_na) {
    assign("encryption", connection_data[5], envir = atsdEnv)
  }
  if (connection_data[6] == "true") {
    cat("Set connection parameters from the", from)
  } else {
    cat("Set connection parameters to NA's")
  }
}
