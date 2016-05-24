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
make_flat <- function(a_list) {
  tags <- a_list$tags
  if (length(names(tags) > 0)) {
    names(tags) <- paste0("tags.", names(tags))
  }
  a_list$tags <- NULL
  data.frame(c(a_list, tags), stringsAsFactors = FALSE)
}
