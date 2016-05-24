##
#  Copyright (c) 2014-2015 LabKey Corporation
#
#  Licensed under the Apache License, Version 2.0 (the "License");
#  you may not use this file except in compliance with the License.
#  You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
##
PACKAGE_ENV = new.env()

# Helper function to get the base curl options used for all http or https requests
#
labkey.setCurlOptions <- function(...)
{
    # default curl options
    options <- curlOptions(ssl.verifyhost=2, ssl.verifypeer=TRUE, followlocation=TRUE, sslversion=1L)

    # check if a certificate bundle has been specified from the environment variable
    vars <- Sys.getenv("RLABKEY_CAINFO_FILE")
    if (nchar(vars[1]) > 0)
    {
        options <- curlOptions(cainfo = vars[1], .opts=c(options))
    }

    # merge in any overrides
    options <- curlOptions(..., .opts=c(options))
    assign("RLABKEY_CURL_OPTIONS", options, envir=PACKAGE_ENV)

    return(get("RLABKEY_CURL_OPTIONS", envir=PACKAGE_ENV))
}

labkey.acceptSelfSignedCerts <- function()
{
    return(labkey.setCurlOptions(ssl.verifyhost=0, ssl.verifypeer=FALSE))
}