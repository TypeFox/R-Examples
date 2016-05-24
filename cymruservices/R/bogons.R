safe_rl <- safely(readLines)

#' Retrieve list of IPv4 "full bogons" from Team Cymru webservice
#'
#' The traditional bogon prefixes (IPV4), plus prefixes that have been allocated to
#' RIRs but not yet assigned by those RIRs to ISPs, end-users, etc. Updated every four hours.
#'
#' Bogons are defined as Martians (private and reserved addresses defined by RFC 1918,
#' RFC 5735, and RFC 6598) and netblocks that have not been allocated to a regional
#' internet registry (RIR) by the Internet Assigned Numbers Authority.
#'
#' Fullbogons are a larger set which also includes IP space that has been allocated to
#' an RIR, but not assigned by that RIR to an actual ISP or other end-user. IANA maintains
#' a convenient IPv4 summary page listing allocated and reserved netblocks, and each RIR
#' maintains a list of all prefixes that they have assigned to end-users. Our bogon reference
#' pages include additional links and resources to assist those who wish to properly filter
#' bogon prefixes within their networks.
#'
#' @param force force a refresh even if the time-frame (4-hours) is not up
#' @param cached_bogons if you pass in the previous result of a call to \code{ipv4_bogoons}
#'        it will be returned if the refresh time constraint has not been met, otherwise
#'        \code{NA} will be returned.
#' @seealso \url{http://www.team-cymru.org/bogon-reference-http.html}
#' @export
#' @examples \dontrun{
#' v4_bogons <- ipv4_bogons()
#' v4_bogons <- ipv4_bogons(cached_bogons=v4_bogons)
#' }
ipv4_bogons <- function(force=FALSE, cached_bogons=NA) {

  BOGONS_URL <- "http://www.team-cymru.org/Services/Bogons/fullbogons-ipv4.txt"

  if (Sys.getenv("CYMRU_LAST_V4_BOGON") == "") {

    bogons <- safe_rl(BOGONS_URL)

    if (is.null(bogons$result)) {
      message("Could not reach team-cymru.org")
      return(NA_character_)
    }

    bogons <- bogons$result

    last_updated <- as.numeric(str_match_all(bogons[1],
                                             "updated ([[:digit:]]+) ")[[1]][,2])

    Sys.setenv(CYMRU_LAST_V4_BOGON=last_updated)

    return(tail(bogons, -1))

  } else {

    delta <- (as.numeric(Sys.time()) - as.numeric(Sys.getenv("CYMRU_LAST_V4_BOGON")))

    if ((delta > 14400) | ((delta < 14400) & force)) {

      bogons <- safe_rl(BOGONS_URL)

      if (is.null(bogons$result)) {
        message("Could not reach team-cymru.org")
        return(NA_character_)
      }

      bogons <- bogons$result

      last_updated <- as.numeric(str_match_all(bogons[1],
                                               "updated ([[:digit:]]+) ")[[1]][,2])

      Sys.setenv(CYMRU_LAST_V4_BOGON=last_updated)

      return(tail(bogons, -1))

    } else {

      message("It has not been 4 hours since the last refresh. Use 'force=TRUE' to force a refresh.")

      return(cached_bogons)

    }

  }

}

#' Retrieve list of IPv6 "full bogons" from Team Cymru webservice
#'
#' IPv6 "fullbogons", all IPv6 prefixes that have not been allocated to RIRs and
#' that have not been assigned by RIRs to ISPs, end-users, etc. Updated every four hours.
#'
#' Bogons are defined as Martians (private and reserved addresses defined by RFC 1918,
#' RFC 5735, and RFC 6598) and netblocks that have not been allocated to a regional
#' internet registry (RIR) by the Internet Assigned Numbers Authority.
#'
#' Fullbogons are a larger set which also includes IP space that has been allocated to
#' an RIR, but not assigned by that RIR to an actual ISP or other end-user. IANA maintains
#' a convenient IPv4 summary page listing allocated and reserved netblocks, and each RIR
#' maintains a list of all prefixes that they have assigned to end-users. Our bogon reference
#' pages include additional links and resources to assist those who wish to properly filter
#' bogon prefixes within their networks.
#'
#' @param force force a refresh even if the time-frame (4-hours) is not up
#' @param cached_bogons if you pass in the previous result of a call to \code{ipv6_bogoons}
#'        it will be returned if the refresh time constraint has not been met, otherwise
#'        \code{NA} will be returned.
#' @seealso \url{http://www.team-cymru.org/bogon-reference-http.html}
#' @export
#' @examples \dontrun{
#' v6_bogons <- ipv6_bogons()
#' v6_bogons <- ipv6_bogons(cached_bogons=v6_bogons)
#' }
ipv6_bogons <- function(force=FALSE, cached_bogons=NA) {

  FULL_BOGONS_URL <- "http://www.team-cymru.org/Services/Bogons/fullbogons-ipv6.txt"

  if (Sys.getenv("CYMRU_LAST_V6_BOGON") == "") {

    bogons <- safe_rl(FULL_BOGONS_URL)

    if (is.null(bogons$result)) {
      message("Could not reach team-cymru.org")
      return(NA_character_)
    }

    bogons <- bogons$result

    last_updated <- as.numeric(str_match_all(bogons[1],
                                             "updated ([[:digit:]]+) ")[[1]][,2])

    Sys.setenv(CYMRU_LAST_V6_BOGON=last_updated)

    return(tail(bogons, -1))

  } else {

    delta <- (as.numeric(Sys.time()) - as.numeric(Sys.getenv("CYMRU_LAST_V6_BOGON")))

    if ((delta > 14400) | ((delta < 14400) & force)) {

      bogons <- safe_rl(FULL_BOGONS_URL)

      if (is.null(bogons$result)) {
        message("Could not reach team-cymru.org")
        return(NA_character_)
      }

      bogons <- bogons$result

      last_updated <- as.numeric(str_match_all(bogons[1],
                                               "updated ([[:digit:]]+) ")[[1]][,2])

      Sys.setenv(CYMRU_LAST_V6_BOGON=last_updated)

      return(tail(bogons, -1))

    } else {

      message("It has not been 4 hours since the last refresh. Use 'force=TRUE' to force a refresh.")

      return(cached_bogons)

    }

  }
}
