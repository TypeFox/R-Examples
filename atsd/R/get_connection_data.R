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
get_connection_data <- function(config_file = NULL) {
  if (is.null(config_file)) {
    config_file <- system.file("connection.config", package = "atsd")
    from <- "the configuration file"
  } else {
    from <- "provided file"
  }
  con <- file(config_file)
  on.exit(close(con))
  lines <- tryCatch(
    {
      readLines(con)
    },
    error=function(er) {
      # message("Here is the original error message:")
      # message(er)
      return(NULL)
    }
  )
  if (config_file == "" || is.null(lines)) {
    message(paste("It seems that ", from, "does not exist."))
    #if (from == "the configuration file") {
    #  message("You could recreate configuration file by means of the save_connection() function.")
    #}
    return(c(NA, NA, NA, NA, NA, "false"))
  }
  url <- get_entry(lines, "url")
  # delete last / if it is
  url <- sub("/$", "", url)
  user <- get_entry(lines, "user")
  password <- get_entry(lines, "password")
  verify <- tolower(get_entry(lines, "verify"))
  if (!is.na(verify)) {
    if (verify %in% c("yes", "true")) {
      verify <- "yes"
    } else {
      verify <- "no"
    }
  }
  encryption <- get_entry(lines, "encryption")
  if (!is.na(encryption)) {
    encryption <- tolower(encryption)
  }
  return(c(url, user, password, verify, encryption, "true"))
}

get_entry <- function(strings, key) {
  line_number <- 1
  value <- NA
  while (line_number <= length(strings) && is.na(value)) {
    line <- gsub("[[:space:]]", "", strings[line_number])
    first_equality <- gregexpr(pattern = "=", line)[[1]][1]
    if (first_equality > 1 && substr(line, 1, first_equality[1] - 1) == key && 
          first_equality < nchar(line)) {
      value <- substr(line, first_equality[1] + 1, nchar(line))
      break
    }
    line_number <- line_number + 1
  }
  return(value)
}
