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
check_time_interval_spelling <- function(ti_name, ti_value) {
  if (!grepl(pattern = "^(-?)[123456789]([[:digit:]]*)[-](Second|Minute|Hour|Day|Week|Month|Quarter|Year)$", 
             x = ti_value, ignore.case = TRUE)) {
    ti_spelling <- paste0(
      "\n",
      "It seems that there is a spelling error in ",
      ti_name,
      " argument of the function.",
      "\n",
      "Provided ",
      ti_name,
      " is:",
      "\n",
      ti_value, 
      "\n",
      'But it should have the form "n-Unit", where n is a number,',
      "\n",
      "and Unit is one of: Second, Minute, Hour, Day, Week, Month, Quarter, Year.",
      "\n",
      'E.g. ',
      ti_name,
      ' = "3-Week", or ',
      ti_name,
      ' = "12-Hour".',
      "\n",
      "Nevertheless try proceed with your query.",
      "\n\n")
  } else {
    ti_spelling <- ""
  }
  return(ti_spelling)
}
