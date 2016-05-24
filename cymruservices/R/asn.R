# handy
sock <- safely(socketConnection)

#' Retrieves BGP Origin ASN info for a list of IPv4 addresses
#'
#' @param ips vector of IPv4 address (character - dotted-decimal)
#' @param timeout	numeric: the timeout (in seconds) to be used for this connection.
#'        Beware that some OSes may treat very large values as zero: however the
#'        POSIX standard requires values up to 31 days to be supported.
#' @return data frame of BGP Origin ASN lookup results
#'   \itemize{
#'     \item \code{as} - AS #
#'     \item \code{ip} - IPv4 (passed in)
#'     \item \code{bgp_refix} - BGP CIDR
#'     \item \code{cc} - Country code
#'     \item \code{registry} - Registry it falls under
#'     \item \code{allocated} - date it was allocated
#'     \item \code{as_ame} - AS name
#'   }
#'   If a socket connection cannot be made (i.e. a network problem on your
#'   end or a service/network problem on their end), all columns will be
#'   \code{NA}.
#' @note The Team Cymru's service is NOT a GeoIP service! Do not use this
#'       function for that as your results will not be accurate.
#'       Data is updated every 4 hours. Also,
#'       A direct connection to TCP Port 43 (WHOIS) is required for most of these
#'       API functions to work properly.
#' @seealso \url{http://www.team-cymru.org/IP-ASN-mapping.html}
#' @export
#' @examples \dontrun{
#' bulk_origin(c("68.22.187.5", "207.229.165.18", "198.6.1.65"))
#' }
bulk_origin <- function(ips, timeout=getOption("timeout")) {

  host <- "v4.whois.cymru.com"
  port <- 43

  # setup query
  cmd <- "begin\nverbose\n"
  ips_c <- paste(unlist(ips), collapse="\n")
  cmd <- sprintf("%s%s\nend\n", cmd, ips_c)

  # setup connection and post query
  con <- sock(host=host, port=port, blocking=TRUE, open="r+", timeout=timeout)
  if (is.null(con$result)) {
    message("Error opening connection to v4.whois.cymru.com")
    return(data.frame(as=rep(NA, length(ips)),
                      ip=rep(NA, length(ips)),
                      bgp_prefix=rep(NA, length(ips)),
                      cc=rep(NA, length(ips)),
                      registry=rep(NA, length(ips)),
                      allocated=rep(NA, length(ips)),
                      as_name=rep(NA, length(ips))))
  }

  con <- con$result
  cat(cmd, file=con)
  response <- readLines(con)
  close(con)

  if (length(response) == 0) {
    message("Error reading from connection to v4.whois.cymru.com")
    return(data.frame(as=rep(NA, length(ips)),
                      ip=rep(NA, length(ips)),
                      bgp_prefix=rep(NA, length(ips)),
                      cc=rep(NA, length(ips)),
                      registry=rep(NA, length(ips)),
                      allocated=rep(NA, length(ips)),
                      as_name=rep(NA, length(ips))))
  }

  # trim header, split fields and convert results
  response <- trim_df(read.csv(textConnection(tail(response, -1)),
                               stringsAsFactors=FALSE, sep="|", header=FALSE))
  names(response) <- c("as", "ip", "bgp_prefix", "cc",
                       "registry", "allocated", "as_name")
  response[response=="NA"] <- NA

  return(response)

}

#' Retrieves BGP Peer ASN info for a list of IPv4 addresses
#'
#' @param ips vector of IPv4 address (character - dotted-decimal)
#' @param timeout	numeric: the timeout (in seconds) to be used for this connection.
#'        Beware that some OSes may treat very large values as zero: however the
#'        POSIX standard requires values up to 31 days to be supported.
#' @return data frame of BGP Peer ASN lookup results
#'   \itemize{
#'     \item \code{peer_as} - peer AS #
#'     \item \code{ip} - IPv4 (passsed in)
#'     \item \code{bgp_prefix} - BGP CIDR block
#'     \item \code{cc} - Country code
#'     \item \code{registry} - Registry it falls under
#'     \item \code{allocated} - date allocated
#'     \item \code{peer_as_name} - peer name
#'   }
#'   If a socket connection cannot be made (i.e. a network problem on your
#'   end or a service/network problem on their end), all columns will be
#'   \code{NA}.
#' @note The Team Cymru's service is NOT a GeoIP service! Do not use this
#'       function for that as your results will not be accurate.
#'       Data is updated every 4 hours. Also,
#'       A direct connection to TCP Port 43 (WHOIS) is required for most of these
#'       API functions to work properly.
#' @seealso \url{http://www.team-cymru.org/IP-ASN-mapping.html}
#' @export
#' @examples \dontrun{
#' bulk_peer(c("68.22.187.5", "207.229.165.18", "198.6.1.65"))
#' }
bulk_peer <- function(ips, timeout=getOption("timeout")) {

  host <- "v4-peer.whois.cymru.com"
  port <- 43

  # setup query
  cmd <- "begin\nverbose\n"
  ips_c <- paste(unlist(ips), collapse="\n")
  cmd <- sprintf("%s%s\nend\n", cmd, ips_c)

  # setup connection and post query
  con <- sock(host=host, port=port, blocking=TRUE, open="r+", timeout=timeout)
  if (is.null(con$result)) {
    message("Error opening connection to v4-peer.whois.cymru.com")
    return(data.frame(peer_as=rep(NA, length(ips)),
                      ip=rep(NA, length(ips)),
                      bgp_prefix=rep(NA, length(ips)),
                      cc=rep(NA, length(ips)),
                      registry=rep(NA, length(ips)),
                      allocated=rep(NA, length(ips)),
                      peer_as_name=rep(NA, length(ips))))
  }

  con <- con$result
  cat(cmd, file=con)
  response <- readLines(con)
  close(con)

  if (length(response) == 0) {
    message("Error reading from connection to v4-peer.whois.cymru.com")
    return(data.frame(peer_as=rep(NA, length(ips)),
                      ip=rep(NA, length(ips)),
                      bgp_prefix=rep(NA, length(ips)),
                      cc=rep(NA, length(ips)),
                      registry=rep(NA, length(ips)),
                      allocated=rep(NA, length(ips)),
                      peer_as_name=rep(NA, length(ips))))
  }


  # trim header, split fields and convert results
  response <- trim_df(read.csv(textConnection(tail(response, -1)),
                               stringsAsFactors=FALSE, sep="|", header=FALSE))
  names(response) <- c("peer_as", "ip", "bgp_prefix", "cc",
                       "registry", "allocated", "peer_as_name")
  response[response=="NA"] <- NA

  return(response)

}

#' Retrieves BGP Origin ASN info for a list of ASN ids
#'
#' @param asns character vector of ASN ids (character)
#' @param timeout	numeric: the timeout (in seconds) to be used for this connection.
#'        Beware that some OSes may treat very large values as zero: however the
#'        POSIX standard requires values up to 31 days to be supported.
#' @return data frame of BGP Origin ASN lookup results
#'   \itemize{
#'     \item \code{as} - AS #
#'     \item \code{cc} - Country code
#'     \item \code{registry} - registry it falls under
#'     \item \code{allocated} - when it was allocated
#'     \item \code{as_name} - name associated with the allocation
#'   }
#'   If a socket connection cannot be made (i.e. a network problem on your
#'   end or a service/network problem on their end), all columns will be
#'   \code{NA}.
#' @note The Team Cymru's service is NOT a GeoIP service! Do not use this
#'       function for that as your results will not be accurate.
#'       Data is updated every 4 hours. Also,
#'       A direct connection to TCP Port 43 (WHOIS) is required for most of these
#'       API functions to work properly.
#' @seealso \url{http://www.team-cymru.org/IP-ASN-mapping.html}
#' @export
#' @examples \dontrun{
#' bulk_origin_asn(c(22822, 1273, 2381, 2603, 2914, 3257, 3356, 11164,
#'                   174, 286, 1299, 2914, 3257, 3356, 3549, 22822))
#' }
bulk_origin_asn <- function(asns, timeout=getOption("timeout")) {

  host <- "v4.whois.cymru.com"
  port <- 43

  # setup query
  cmd <- "begin\nverbose\n"
  ips <- paste(unlist(ifelse(grepl("^AS", asns), asns,
                             sprintf("AS%s", asns))), collapse="\n")
  cmd <- sprintf("%s%s\nend\n", cmd, ips)

  # setup connection and post query
  con <- sock(host=host, port=port, blocking=TRUE, open="r+", timeout=timeout)
  if (is.null(con$result)) {
    message("Error opening connection to v4.whois.cymru.com")
    return(data.frame(as=rep(NA, length(asns)),
                      cc=rep(NA, length(asns)),
                      registry=rep(NA, length(asns)),
                      allocated=rep(NA, length(asns)),
                      as_name=rep(NA, length(asns))))
  }

  con <- con$result
  cat(cmd, file=con)
  response <- readLines(con)
  close(con)

  if (length(response) == 0) {
    message("Error reading from connection to v4.whois.cymru.com")
    return(data.frame(as=rep(NA, length(asns)),
                      cc=rep(NA, length(asns)),
                      registry=rep(NA, length(asns)),
                      allocated=rep(NA, length(asns)),
                      as_name=rep(NA, length(asns))))
  }

  # trim header, split fields and convert results
  response <- trim_df(read.csv(textConnection(tail(response, -1)),
                               stringsAsFactors=FALSE, sep="|", header=FALSE))
  names(response) <- c("as", "cc", "registry", "allocated", "as_name")
  response[response=="NA"] <- NA

  return(response)

}
