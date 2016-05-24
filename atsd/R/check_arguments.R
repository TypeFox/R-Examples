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
check_arguments <- function(export_type, metric, selection_interval) {
  if (!(tolower(export_type) %in% c("history", "forecast"))) {
    et_mistake <- paste0(
      "\n",
      'Export_type argument should be one of: "History" or "Forecast"!',
      "\n")
  } else {
    et_mistake <- ""
  }
  if (missing(metric)) {
    m_missing <- paste0(
      "\n",
      "Your should provide metric argument to the function.",
      "\n",
      "To view available metrics use the 'get_metrics' function.",
      "\n")
  } else {
    m_missing <- ""
  }
  if (missing(selection_interval)) {
    si_missing <- paste0(
      "\n",
      "Your should provide selection_interval argument to the function.",
      "\n",
      'It should have the form "n-Unit", where n is a number,',
      "\n",
      "and Unit is one of: Second, Minute, Hour, Day, Week, Month, Quarter, Year.",
      "\n",
      'E.g. selection_interval = "3-Week", or selection_interval = "12-Hour".',
      "\n\n")
  } else {
    si_missing <- ""
  }
  return(paste0(et_mistake, m_missing, si_missing))
}
