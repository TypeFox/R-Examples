
getFlights <- function(status = c("ankunft", "abflug"),
    hour = formatC(c(0, 3:11 * 2), flag = "0", width = 2, format = "d")) {

  stopifnot(require("XML"))

  status <- match.arg(status)
  hour <- match.arg(hour, several.ok = TRUE)

  URLs <- switch(status,
      "ankunft" =  paste("http://www.munich-airport.de/de/consumer/fluginfo/ankunft/h",
                         hour, "00_de_L.jsp", sep = ""),
      "abflug" = paste("http://www.munich-airport.de/de/consumer/fluginfo/abflug/h",
                       hour, "00_de_S.jsp", sep = ""))

  websites <- lapply(URLs, htmlTreeParse, useInternalNodes = TRUE, encoding = "UTF-8")

  ## raw Data List
  rDL<- lapply(websites, xpathSApply,
               path="//table//input[@type='hidden']", xmlAttrs)

  lapply(websites, free)

  ### remove pages without a flight table
  rDL <- rDL[sapply(rDL, length) > 0]

  dataList <- lapply(rDL, function(z) {
                  ret <- matrix(z["value", ], ncol = 20, byrow = TRUE)
                  colnames(ret) <- z["name", 1:20]
                  as.data.frame(ret, stringsAsFactors = FALSE)
              })

  ## ornamente
  resDF <- do.call(rbind, dataList)
  resDF <- resDF[, !colnames(resDF) %in% c("key", "zurueck")]
  resDF$ett[resDF$ett == ""] <- resDF$stt[resDF$ett == ""]
  ## resDF$stt <- strftime(resDF$stt, format = "%Y-%m-%d %H:%M:%S.0")
  ## resDF$ett <- strftime(resDF$ett, format = "%Y-%m-%d %H:%M:%S.0")
  return(resDF)
}

# muc <- rbind(getFlights("ankunft"), getFlights("abflug"))
# save(muc, file = paste("../data/", Sys.Date(), "_Flugdaten", ".RData", sep = "" ))
