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
set_https_options <- function() {
  if (substr(get("url", envir = atsdEnv), 1, 5) == "https") {
    encryption_code <- switch(get("encryption", envir = atsdEnv),
                              "tls1" = RCurl::SSLVERSION_TLSv1,
                              "ssl2" = RCurl::SSLVERSION_SSLv2,
                              "ssl3" = RCurl::SSLVERSION_SSLv3,
                              RCurl::SSLVERSION_DEFAULT
    )
    return(list(ssl.verifypeer = get("verify", envir = atsdEnv), sslversion = encryption_code))
  }
  return(list())
}
