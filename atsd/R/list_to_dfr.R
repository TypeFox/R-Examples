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
list_to_dfr <- function(lst) {
  col_names <- character()
  invisible(lapply(lst, function(x) {col_names <<- union(col_names, names(x))}))
  df <- data.frame(matrix(NA, nrow = length(lst), ncol = length(col_names)))
  colnames(df) <- col_names
  counter <- 1
  f <- function(x) {
    for (field in names(x)) {
      df[counter ,field] <<- x[[field]]
    }
    counter <<- counter + 1
  }
  invisible(lapply(lst, f))
  return(df)
}