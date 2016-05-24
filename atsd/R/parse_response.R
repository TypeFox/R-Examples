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
parse_response <- function(body, header) {
  is_OK <- (gregexpr(pattern = "200 OK", header)[[1]] != -1)  
  if (!is_OK) {
    message <- paste0(
    "\n",
    "There is some problem with processing of your request.",
    "\n",
    "The first line of ATSD server response is:",
    "\n")
    line_ends <- gregexpr(pattern = "\n", header)[[1]]
    if (line_ends[1] == -1) {
      message <- paste0(message, header, "\n\n")
    } else {
      message <- paste0(message, substr(header, 1, line_ends[1] - 1), "\n\n")
    }
    dfr <- data.frame()
  } else {
    message <- paste0("\nYour request has succeeded. ATSD server response is: 200 OK.")
    if (nchar(body) == 0) {
      message <- paste0(message,
        "\n",
        "ATSD could not find any data for your query and returns empty data set!",
        "\n",
        "For the help with arguments of the query look at the function help page.")
        dfr <- data.frame()
    } else {
      t_con <- textConnection(body)
      rm(body)
      dfr <- read.csv(t_con, stringsAsFactors = FALSE)
      if (nrow(dfr) > 0)
      dfr$Timestamp <- strptime(as.character(dfr$Timestamp), format='%Y-%m-%d %H:%M:%S', tz = 'GMT')
      close(t_con)
      #       message <- paste0(message,
      #                         "\n",
      #                         "Head of the fetched data frame. \n\n",
      #                         head(dfr),
      #                         "\n\n",
      #                         "Summary of the data frame. \n\n",
      #                         summary(dfr))
    }
  }
  return(list(dfr, message))
}
  