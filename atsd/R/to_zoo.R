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
#' Build zoo object from data frame.
#' 
#' @description
#' The function builds a zoo object based on the 'Timestamp' and 'Value' 
#' columns of the given data frame. Information from other columns will be lost.
#' To use this function the 'zoo' package should be installed.
#' To install the 'zoo' package type:
#' \code{install.packages("zoo")}.
#' 
#' @param dfr 
#' The data frame with 'Timestamp' and 'Value' columns.
#' 
#' @export
to_zoo <- function(dfr) {
  if (!requireNamespace("zoo", quietly = TRUE)) {
    msg <- '"zoo"package needed for this function to work. \n'
    msg <- paste0(msg, 'Please install it with: install.packages("zoo")')
    stop(msg, call. = FALSE)
  }
  if (all(c("Timestamp", "Value") %in% names(dfr))) {
    return(zoo::zoo(dfr$Value, dfr$Timestamp))
  } else {
    cat("Can not convert given data frame to zoo object.")
    cat("\nThe data frame should has 'Timestamp' and 'Value' columns.")
  }
}
