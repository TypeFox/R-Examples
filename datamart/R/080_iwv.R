#' website traffic from IWV
#'
#' website traffic as tracked by iwv online.
#' Use yyyymm for query() as resource.
#' This is an example for the \code{urldata} function.
#'
#' @seealso \code{\link{urldata}}
#'
#' @return UrlData object
#' @references 
#' \href{http://en.wikipedia.org/wiki/Informationsgemeinschaft_zur_Feststellung_der_Verbreitung_von_Werbetraegern}{Wikipedia}
#' @docType data
#' @export
iwv_online <- function() urldata(
  resource="Iwv",
  template="http://ausweisung.ivw-online.de/i.php?s=$(detail_level)&mz=$(year)$(month)&csv=1", # s defines the detail level 1:3
  detail_level=1,
  year=strftime(Sys.time(), "%Y"),
  month=strftime(Sys.time(), "%m"), 
  extract.fct=readLines,
  transform.fct=function(x) {
    i <- min(which(x==""))
    dat <- read.csv2(textConnection(x[1:(i-1)]), header=TRUE, na.strings="---", stringsAsFactors=FALSE)
    dat$X <- NULL
    dat$Hinweis <- as.character(dat$Hinweis)
    dat[is.na(dat$Hinweis), "Hinweis"] <- ""

    numCols <- c("Visits.Gesamt", "Visits.Inland", "Visits.Ausland", "Visits.Vormonat", "Visits.Vorjahresmonat",  "Redaktioneller.Content.Kategorienvisits.Gesamt", "Redaktioneller.Content.Kategorienvisits.Inland", "Redaktioneller.Content.Kategorienvisits.Ausland", "User.generierter.Content.Kategorienvisits.Gesamt", "User.generierter.Content.Kategorienvisits.Inland", "User.generierter.Content.Kategorienvisits.Ausland", "E.Commerce.Kategorienvisits.Gesamt", "E.Commerce.Kategorienvisits.Inland", "E.Commerce.Kategorienvisits.Ausland", "Kommunikation.Kategorienvisits.Gesamt", "Kommunikation.Kategorienvisits.Inland", "Kommunikation.Kategorienvisits.Ausland", "Suchmaschinen..Verzeichnisse.und.Auskunftsdienste.Kategorienvisits.Gesamt", "Suchmaschinen..Verzeichnisse.und.Auskunftsdienste.Kategorienvisits.Inland", "Suchmaschinen..Verzeichnisse.und.Auskunftsdienste.Kategorienvisits.Ausland", "Spiele.Kategorienvisits.Gesamt", "Spiele.Kategorienvisits.Inland", "Spiele.Kategorienvisits.Ausland", "Diverses.Kategorienvisits.Gesamt", "Diverses.Kategorienvisits.Inland", "Diverses.Kategorienvisits.Ausland")
    for (clmn in numCols) dat[,clmn] <- as.numeric(gsub("\\.", "", dat[,clmn]))
    
    prozCols <- colnames(dat)[grep("Proz", colnames(dat))]
    for (clmn in prozCols) {
      tmp <- gsub(" %", "", dat[,clmn])
      tmp <- gsub(",", ".", tmp)
      dat[,clmn] <- as.numeric(tmp)
    }

    return(dat)
  }
)

