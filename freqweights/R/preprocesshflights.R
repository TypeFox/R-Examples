## Emilio Torres Manzanera
## University of Oviedo
## Time-stamp: <2015-06-21 09:20 emilio on emilio-despacho>
## ============================================================


##' Preprocess the \code{hflights} data.
##' 
##' Preprocess the \code{hflights} data to get a nice format. 
##' 
##' Preprocess this data set.
##' 
##' @param hflights \code{hflights} data set 
##' @return A preprocessed data set. 
##' @seealso \code{\link[hflights]{hflights}} 
##' @keywords manip
##' @export
##' @import dplyr
##' @importFrom data.table fread setnames
##' @examples
##' 
##' if(require(hflights)) {
##' a <- preprocesshflights(hflights[1:10000,])
##' summary(a)
##' }
##' \dontrun{
##' library(hflights)
##' ## We create a file with no header
##' input <- "hflights.csv"
##' write.table(hflights,file=input,sep=",",
##'             row.names=FALSE,col.names=FALSE)
##' ## Inefficient way to read the data. Just as example
##' lines <- readLines(input)
##' lines <- gsub("\"","",lines,fixed=TRUE )
##' x <- strsplit(lines,",")
##' x <- as.data.frame(do.call("rbind",x)) 
##' x <- preprocesshflights(x)
##' summary(x)
##' }
##' 
preprocesshflights <- function(hflights){
  nams <- c("Year", "Month", "DayofMonth", "DayOfWeek", "DepTime", "ArrTime",
            "UniqueCarrier", "FlightNum", "TailNum", "ActualElapsedTime",
            "AirTime", "ArrDelay", "DepDelay", "Origin", "Dest", "Distance",
            "TaxiIn", "TaxiOut", "Cancelled", "CancellationCode", "Diverted")
  varsfactor <- c("UniqueCarrier","Origin","Dest","CancellationCode")
  listoflevels <- list()
  listoflevels[["CancellationCode"]] <- c("A", "B", "C", "D")
  listoflevels[["UniqueCarrier"]] <- c("AA", "AS", "B6", "CO", "DL",
                                       "EV", "F9", "FL", "MQ", "OO",
                                       "UA", "US", "WN", "XE", "YV")
  listoflevels[["Origin"]] <-  c("HOU", "IAH")
  listoflevels[["Dest"]] <-
    c("ABQ", "AEX", "AGS", "AMA", "ANC", "ASE", "ATL", "AUS", "AVL",
      "BFL", "BHM", "BKG", "BNA", "BOS", "BPT", "BRO", "BTR", "BWI",
      "CAE", "CHS", "CID", "CLE", "CLT", "CMH", "COS", "CRP", "CRW",
      "CVG", "DAL", "DAY", "DCA", "DEN", "DFW", "DSM", "DTW", "ECP",
      "EGE", "ELP", "EWR", "FLL", "GJT", "GPT", "GRK", "GRR", "GSO",
      "GSP", "GUC", "HDN", "HNL", "HOB", "HRL", "HSV", "IAD", "ICT",
      "IND", "JAN", "JAX", "JFK", "LAS", "LAX", "LBB", "LCH", "LEX",
      "LFT", "LGA", "LIT", "LRD", "MAF", "MCI", "MCO", "MDW", "MEM",
      "MFE", "MIA", "MKE", "MLU", "MOB", "MSP", "MSY", "MTJ", "OAK",
      "OKC", "OMA", "ONT", "ORD", "ORF", "PBI", "PDX", "PHL", "PHX",
      "PIT", "PNS", "PSP", "RDU", "RIC", "RNO", "RSW", "SAN", "SAT",
      "SAV", "SDF", "SEA", "SFO", "SHV", "SJC", "SJU", "SLC", "SMF",
      "SNA", "STL", "TPA", "TUL", "TUS", "TYS", "VPS", "XNA")
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  ftime <- function(x){
    x <- formatC(x, width = 4, format = "d", flag = "0")
    x[ x=="  NA"] <- NA
    strftime(strptime(x,"%H%M"),"%H:%M")
  }
  ffactor <- function(x, levels){
    factor(as.character(x),levels=levels)
  }
  ### To avodi R CMD CHECK warnings
  DepTime <-  ArrTime <- UniqueCarrier <- Origin <- Dest <- CancellationCode <- Year <- NULL
  Year <-  Month <- DayOfWeek <- DayofMonth <- ActualElapsedTime <- AirTime <- NULL
  ArrDelay <-  DepDelay <- Distance <- TaxiOut <- TaxiIn <- NULL
#####
  setnames(hflights,nams)
  ##print(colnames(hflights))
  hflights %>% mutate(DepTime = ftime(DepTime),
                      ArrTime = ftime(ArrTime)) %>%
                        mutate(UniqueCarrier = ffactor( UniqueCarrier, listoflevels[["UniqueCarrier"]]),
                               Origin = ffactor( Origin, listoflevels[["Origin"]]),
                               Dest = ffactor( Dest, listoflevels[["Dest"]]),
                               CancellationCode = ffactor( CancellationCode, listoflevels[["CancellationCode"]])
                               ) %>%
                                 mutate(Year = as.integer(as.character(Year)),
                                        Month = as.integer(as.character(Month)),
                                        DayofMonth = as.integer(as.character(DayofMonth)),
                                        DayOfWeek = as.integer(as.character(DayOfWeek)),
                                        ActualElapsedTime = as.integer(as.character(ActualElapsedTime)),
                                        AirTime =  as.integer(as.character(AirTime)),
                                        ArrDelay =  as.integer(as.character(ArrDelay)),
                                        DepDelay =  as.integer(as.character(DepDelay)),
                                        Distance =  as.integer(as.character(Distance)),
                                        TaxiIn =  as.integer(as.character(TaxiIn)),
                                        TaxiOut =  as.integer(as.character(TaxiOut))
                                        )
}
