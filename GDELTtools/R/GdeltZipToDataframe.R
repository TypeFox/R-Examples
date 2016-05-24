GdeltZipToDataframe <- function(f, daily=FALSE, verbose=TRUE) {
  
  gdelt.colNames <- c("GLOBALEVENTID", "SQLDATE", "MonthYear", "Year", "FractionDate", "Actor1Code", "Actor1Name", "Actor1CountryCode", "Actor1KnownGroupCode", "Actor1EthnicCode", "Actor1Religion1Code", "Actor1Religion2Code", "Actor1Type1Code", "Actor1Type2Code", "Actor1Type3Code", "Actor2Code", "Actor2Name", "Actor2CountryCode", "Actor2KnownGroupCode", "Actor2EthnicCode", "Actor2Religion1Code", "Actor2Religion2Code", "Actor2Type1Code", "Actor2Type2Code", "Actor2Type3Code", "IsRootEvent", "EventCode", "EventBaseCode", "EventRootCode", "QuadClass", "GoldsteinScale", "NumMentions", "NumSources", "NumArticles", "AvgTone", "Actor1Geo_Type", "Actor1Geo_FullName", "Actor1Geo_CountryCode", "Actor1Geo_ADM1Code", "Actor1Geo_Lat", "Actor1Geo_Long", "Actor1Geo_FeatureID", "Actor2Geo_Type", "Actor2Geo_FullName", "Actor2Geo_CountryCode", "Actor2Geo_ADM1Code", "Actor2Geo_Lat", "Actor2Geo_Long", "Actor2Geo_FeatureID", "ActionGeo_Type", "ActionGeo_FullName", "ActionGeo_CountryCode", "ActionGeo_ADM1Code", "ActionGeo_Lat", "ActionGeo_Long", "ActionGeo_FeatureID", "DATEADDED")
  gdelt.colClasses <- c(GLOBALEVENTID="integer", SQLDATE="integer", MonthYear="integer", Year="integer", FractionDate="numeric", Actor1Code="character", Actor1Name="character", Actor1CountryCode="character", Actor1KnownGroupCode="character", Actor1EthnicCode="character", Actor1Religion1Code="character", Actor1Religion2Code="character", Actor1Type1Code="character", Actor1Type2Code="character", Actor1Type3Code="character", Actor2Code="character", Actor2Name="character", Actor2CountryCode="character", Actor2KnownGroupCode="character", Actor2EthnicCode="character", Actor2Religion1Code="character", Actor2Religion2Code="character", Actor2Type1Code="character", Actor2Type2Code="character", Actor2Type3Code="character", IsRootEvent="integer", EventCode="character", EventBaseCode="character", EventRootCode="character", QuadClass="integer", GoldsteinScale="numeric", NumMentions="integer", NumSources="integer", NumArticles="integer", AvgTone="numeric", Actor1Geo_Type="integer", Actor1Geo_FullName="character", Actor1Geo_CountryCode="character", Actor1Geo_ADM1Code="character", Actor1Geo_Lat="numeric", Actor1Geo_Long="numeric", Actor1Geo_FeatureID="integer", Actor2Geo_Type="integer", Actor2Geo_FullName="character", Actor2Geo_CountryCode="character", Actor2Geo_ADM1Code="character", Actor2Geo_Lat="numeric", Actor2Geo_Long="numeric", Actor2Geo_FeatureID="integer", ActionGeo_Type="integer", ActionGeo_FullName="character", ActionGeo_CountryCode="character", ActionGeo_ADM1Code="character", ActionGeo_Lat="numeric", ActionGeo_Long="numeric", ActionGeo_FeatureID="integer", DATEADDED="integer")
  
  daily.cols.w.errors <- c("Actor1Geo_FeatureID", "Actor2Geo_FeatureID", "ActionGeo_FeatureID", "IsRootEvent", "GoldsteinScale", "NumMentions", "AvgTone", "NumSources", "NumArticles", "Actor1Geo_Lat", "Actor1Geo_Long", "Actor2Geo_Lat", "Actor2Geo_Long", "ActionGeo_Lat", "ActionGeo_Long")
  if(daily) {
    gdelt.colNames <- c(gdelt.colNames, "SOURCEURL")
    gdelt.colClasses[daily.cols.w.errors]<-"character"
  }
  if(verbose) cat("Ingesting", f, "\n")
  out <- read.delim(unz(f, unzip(f, list=TRUE)$Name[1]), 
                    header=FALSE, 
                    stringsAsFactors=FALSE,
                    col.names=gdelt.colNames,
                    colClasses=gdelt.colClasses,
                    na.strings=c("NA", ""),
                    strip.white=TRUE)
  
  if(daily) {
    for(field.name in daily.cols.w.errors) {
      suppressWarnings(out[,field.name] <- as.numeric(out[,field.name]))
    }
  } else {
    out <- data.frame(out, "", stringsAsFactors=FALSE)
    names(out) <- c(gdelt.colNames, "SOURCEURL")
  }
  return(out)
}