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
get_request <- function(  export_type,
                          metric,
                          entity = NA,
                          entity_group = NA,
                          tags = character(),  
                          selection_interval,
                          end_time = NA,
                          aggregate_interval = NA,
                          interpolation = "None",
                          aggregate_statistics = character()) {
  request <- '{'
  
  tag_spelling <- ""
  if (length(tags) > 0) {
    request <- paste0(request, '"tags":[')  
    for (tag in tags) {
      equality_entries <- gregexpr(pattern = "=", tag)[[1]]
      if (equality_entries[1] == -1) {
        tag_spelling <- paste0(tag_spelling,
        "\n",
        "There is an error in the tag argument of the function.",
        "\n",
        "One of tag arguments is:",
        "\n",
        tag, 
        "\n",
        'But a tag should have the form "key=value"!',
        "\n",
        "We try proceed with your query without this tag.",
        "\n\n")
        next
      }
      key <- substr(tag, 1, equality_entries[1] - 1)
      value <- substr(tag, equality_entries[1] + 1, nchar(tag))
      tag_string <- paste0('{"k":"', key, '","v":"', value, '"},')
      request <- paste0(request, tag_string)
    }
    request <- paste0(substr(request, 1, nchar(request) - 1), '],')
  }
  
  request <- paste0(request, '"m":"', metric, '",')
  
  if (!is.na(entity)) {
    request <- paste0(request, '"e":"', entity, '",')
  }
  
  if (!is.na(entity_group)) {
    request <- paste0(request, '"g":"', entity_group, '",')
  }
  
  si_spelling <- check_time_interval_spelling("selection_interval", selection_interval)

  request <- paste0(request, '"si":"', toupper(selection_interval), '",')
  
  if (!is.na(end_time)) {
    request <- paste0(request, '"et":"', end_time, '",')
  }
  
  request <- paste0(request, '"t":"', toupper(export_type), '",')
  
  request <- paste0(request, '"f":"CSV",')
  
  ai_missed <- ""
  ai_spelling <- ""
  if (!is.na(aggregate_interval) || length(aggregate_statistics) > 0) {
    if (xor(is.na(aggregate_interval), length(aggregate_statistics) == 0)) {
      ai_missed <- paste0(
        "\n",
        "One of aggregate_interval or aggregate_statistics argument not provided to the function.",
        "\n",
        "So data will be fetched without aggregetion.",
        "\n",
        "If you need aggregation please specify both of aggregate_interval and",
        "\n",
        "aggregate_statistics arguments.",
        "\n\n")
    } else {
      ai_spelling <- check_time_interval_spelling("aggregate_interval", aggregate_interval)
      request <- paste0(request, '"ai":"', toupper(aggregate_interval), '",')
      
      request <- paste0(request, '"a":["')  
      arg_names <- c("Avg", "Min", "Max", "Sum", "Count", "StDev", "WAvg", "WTAvg", "Percentile 50", 
                     "Percentile 75", "Percentile 90", "Percentile 95", "Percentile 99", 
                     "Percentile 99.5", "Percentile 99.9")
      req_names <- c("AVG", "MIN", "MAX", "SUM", "COUNT", "STANDARD_DEVIATION", "WAVG", "WTAVG", 
                     "PERCENTILE_50", "PERCENTILE_75", "PERCENTILE_90", "PERCENTILE_95", 
                     "PERCENTILE_99", "PERCENTILE_995", "PERCENTILE_999")
      for (stat in aggregate_statistics) {
        request <- paste0(request, req_names[which(stat == arg_names)], '","')
      }
      request <- paste0(substr(request, 1, nchar(request) - 2), '],')
      
      if (interpolation != "None") {
        request <- paste0(request, '"i":"', toupper(interpolation), '",')
      }
    }
  }
  request <- paste0(substr(request, 1, nchar(request) - 1), '}')
  warnings <- paste0(tag_spelling, si_spelling, ai_missed, ai_spelling)
  
  return(c(request, warnings))
}
