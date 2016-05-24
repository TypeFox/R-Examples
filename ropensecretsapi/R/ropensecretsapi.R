#' @import RJSONIO RCurl
#'
NULL

rOpenSecretsAPI.env <- new.env()

#' Allows the user to set the OpenSecrets.org API key once thereby removing the
#' need to pass it in for each function call in this package.
#'
#' @param apiKey The OpenSecrets.org API key.
#'
#' @examples
#'  SetAPIKey ("Example API Key")
#'
#' @export
#'
SetAPIKey <- function (apiKey) {
    assign("openSecretsAPIKey", apiKey, envir = rOpenSecretsAPI.env)
}

#' Function is used to add the key and value to the params that will be sent on
#' every call to OpenSecrets.org; an exception will be thrown if either of
#' these parameters is null.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#' @param key The key to add to the params variable.
#' @param value The value to add to the params variable.
#'
.AddKeyValuePair <- function (params, key, value) {

    if (is.null (params)) {
        stop ("The params parameter cannot be null.")
    }

    if (is.null (key)) {
        stop ("The key parameter cannot be null.")
    }

    if (is.null (value)) {
        stop ("The value parameter cannot be null.")
    }

    params[key] = value

    return (params)
}

#' Function is used to set the output parameter to JSON.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#'
.AddOutputAsJSON <- function (params) {
    params <- .AddKeyValuePair (params, "output", "json")
    return (params)
}

#' Function is used to set the OpenSecrets.org API method being called.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#'
#' @param method The web method being invoked -- see
#' \href{https://www.opensecrets.org/resources/create/api_doc.php}{here}.
#'
.AddMethod <- function (params, method) {
    params <- .AddKeyValuePair (params, "method", method)
    return (params)
}

#' Function is used to add the API key to the params that will be sent on
#' every call to OpenSecrets.org; an exception will be thrown if this value is
#' null.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#'
.AddAPIKey <- function (params) {
    params <- .AddKeyValuePair (params, "apikey", rOpenSecretsAPI.env$openSecretsAPIKey)
    return (params)
}

#' Function performs the request for data from the OpenSecrets.org website and
#' returns the result as a data frame.
#'
#' @param method The web method being invoked -- see
#' \href{https://www.opensecrets.org/resources/create/api_doc.php}{here}.
#'
#' @param params Any parameter pertaining to a specific web method invocation.
#'
.DoGet <- function (method, params) {

    params <- .AddAPIKey (params)

    params <- .AddOutputAsJSON (params)

    params <- .AddMethod (params, method)

    data <- getForm("http://www.opensecrets.org/api/", .params=params)
    dataFrame <- RJSONIO::fromJSON(data)

    return (dataFrame)
}

#' Provides the top organizations contributing to specified politician.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=candContrib}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (cid="N00007360", cycle="2012")
#'  tryCatch(
#'      candContribData <- GetCandContribData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetCandContribData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=candContrib&cid=N00007360&cycle=2012&apikey=__apikey__
    #
    results <- .DoGet ("candContrib", params)
    return (results)
}

#' Provides total contributed to specified candidate from specified industry
#' for specified cycle.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=candIndByInd}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (cid="N00007360", cycle="2012", ind="K02")
#'  tryCatch(
#'      candIndByIndData <- GetCandIndByIndData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetCandIndByIndData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=CandIndByInd&cid=N00007360&cycle=2012&ind=K02&apikey=_API_KEY_
    #
    results <- .DoGet ("CandIndByInd", params)
    return (results)
}

#' Provides the top industries contributing to a specified politician.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=candIndustry}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (cid="N00007360", cycle="2012", ind="K02")
#'  tryCatch(
#'      candIndustryData <- GetCandIndustryData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetCandIndustryData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=candIndustry&cid=N00007360&cycle=2012&apikey=__apikey__
    #
    results <- .DoGet ("candIndustry", params)
    return (results)
}

#' Provides the top industries contributing to a specified politician.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=candSector}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (cid="N00007360", cycle="2012")
#'  tryCatch(
#'      candSectorData <- GetCandSectorData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#' @export
#'
GetCandSectorData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=candSector&cid=N00007360&cycle=2012&apikey=__apikey__
    #
    results <- .DoGet ("candSector", params)
    return (results)
}

#' Provides summary fundraising information for a specified politician.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=candSummary}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (cid="N00007360", cycle="2012")
#'  tryCatch(
#'      candSummaryData <- GetCandSummaryData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetCandSummaryData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=candSummary&cid=N00007360&cycle=2012&apikey=__apikey__
    #
    results <- .DoGet ("candSummary", params)
    return (results)
}

#' Provides summary fundraising information for a specific committee, industry
#' and congress number.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=congCmteIndus}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (congno="112", indus="F10", cmte="HARM")
#'  tryCatch(
#'      congCmteIndusData <- GetCongCmteIndusData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetCongCmteIndusData <- function (params) {
    #
    # http://www.opensecrets.org/api/?method=congCmteIndus&congno=112&indus=F10&cmte=HARM&apikey=__apikey__
    #
    results <- .DoGet ("congCmteIndus", params)
    return (results)
}

#' Provides a list of legislators and associated attributes.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=getLegislators}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (id="NJ")
#'  tryCatch(
#'      legislatorsData <- GetLegislatorsData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetLegislatorsData <- function (params) {
  #
  # http://www.opensecrets.org/api/?method=getLegislators&id=NJ&apikey=__apikey__
  #
  results <- .DoGet ("getLegislators", params)
  return (results)
}

#' Provides organization data.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=getOrgs}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (org="Microsoft")
#'  tryCatch(
#'      orgsData <- GetOrgsData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetOrgsData <- function (params) {
  #
  # N/A
  #
  results <- .DoGet ("getOrgs", params)
  return (results)
}

#' Returns PFD information for a member of Congress.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=memPFDprofile}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (year="2010", cid="N00007360")
#'  tryCatch(
#'      memPFDprofileData <- GetMemPFDprofileData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#'  }
#'
#' @export
#'
GetMemPFDprofileData <- function (params) {
  #
  # http://www.opensecrets.org/api/?method=memPFDprofile&year=2010&cid=N00007360&output=xml&apikey=__apikey__
  #
  results <- .DoGet ("memPFDprofile", params)
  return (results)
}

#' Provides organization summary data.
#'
#' @param params Any parameter accepted by this web service call -- see
#' \href{https://www.opensecrets.org/api/?output=doc&method=orgSummary}{here}.
#'
#' @examples
#'  \dontrun{
#'  SetAPIKey ("ENTER YOUR PRIVATE API KEY HERE.")
#'  params <- list (id="123")
#'  tryCatch(
#'      orgSummaryData <- GetOrgSummaryData (params),
#'      error =
#'          function (e) {    
#'              print (
#'                  paste (
#'                      "An exception was thrown -- details follow: ",
#'                      e,
#'                      sep=""
#'                  )
#'              )
#'          }
#'      )
#' }
#'
#' @export
#'
GetOrgSummaryData <- function (params) {
  #
  # N/A
  #
  results <- .DoGet ("orgSummary", params)
  return (results)
}