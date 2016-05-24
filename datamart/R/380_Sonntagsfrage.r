#' Election polls for Germany
#'
#' This function exposes an interface to the
#' IfD/Allensbach website's poll data.
#'
#' It serves as a proof of concept for the \code{urldata} function.
#'
#' @param resource       the resource name (character)
#'
#' @seealso \code{\link{urldata}}
#' @export
allensbach <- function(resource="allensbach") urldata(
    resource=resource,
    template="http://www.ifd-allensbach.de/studien-und-berichte/sonntagsfrage/gesamt.html",
    extract.fct=XML::htmlParse,
    transform.fct=function(x) {
        nodes <- try(XML::getNodeSet(x, '//tr[@class="odd" or @class="even"]'), silent=TRUE)
        
        if(inherits(nodes, "try-error")) return(data.frame())
        if(length(nodes)==0) return(data.frame())
        dat <- as.data.frame(t(sapply(nodes, function(n) c(
            timestamp=XML::xmlValue(n[[1]]),
            afd=XML::xmlValue(n[[3]]),
            cducsu=XML::xmlValue(n[[4]]),
            spd=XML::xmlValue(n[[5]]),
            fdp=XML::xmlValue(n[[6]]),
            gruene=XML::xmlValue(n[[7]]),
            linke=XML::xmlValue(n[[8]]),
            piraten=XML::xmlValue(n[[9]]),
            sonstige=XML::xmlValue(n[[10]])
        ))), stringsAsFactors=FALSE)
        
        parse_number <- function(x) {
            idx <- x=="--"
            x[!idx] <- gsub(",", ".", x[!idx], fixed=TRUE)
            x[idx] <- NA
            return(as.numeric(x))
        }
        
        dat$timestamp=strptime(strtail(dat$timestamp, 10), "%d.%m.%Y")
        dat$afd=parse_number(dat$afd)
        dat$cducsu=parse_number(dat$cducsu)
        dat$spd=parse_number(dat$spd)
        dat$fdp=parse_number(dat$fdp)
        dat$gruene=parse_number(dat$gruene)
        dat$linke=parse_number(dat$linke)
        dat$piraten=parse_number(dat$piraten)
        dat$sonstige=parse_number(dat$sonstige)
        
        return(dat)
    }
)

