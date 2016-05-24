# Purpose        : Automated generation of (spatial) metadata
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : Michael Blaschek (blaschek@geographie.uni-kiel.de) and Eloi Ribeiro (eloi.carvalhoribeiro@wur.nl); 
# Dev Status     : Beta
# Note           : Based on the US gov sp metadata standards [http://www.fgdc.gov/metadata/csdgm/], which can be converted to "ISO 19139" XML schema; and on the INSPIRE Metadata Implementing Rules [http://inspire.jrc.ec.europa.eu/documents/Metadata/MD_IR_and_ISO_20131029.pdf]

## internal methods:
setMethod("GetNames", "SpatialMetadata", function(obj){paste(obj@field.names)})
setMethod("GetPalette", "SpatialMetadata", function(obj){obj@palette})

## Generate a spMetadata class object:
.spMetadata.Spatial <- function(
  obj,   
  xml.file, ## optional input metadata file
  out.xml.file,
  md.type = c("FGDC", "INSPIRE")[1],
  generate.missing = TRUE,
  GoogleGeocode = FALSE,
  signif.digit = 3,
  colour_scale,
  color = NULL,
  bounds,
  legend_names,
  icons,
  validate.schema = FALSE,
  ...  
)
{
  
  ## Check if there is a data slot:
  if(!("data" %in% slotNames(obj))){
    stop("Object of class 'Spatial*DataFrame' required.")
  }  

  ## Use metadata file if it does exit:
  if(!missing(xml.file)){
    if(file.exists(xml.file)){
      message(paste("Reading the metadata file: ", xml.file, sep=""))
      ret <- xmlTreeParse(xml.file)  
      ml <- xmlRoot(ret)
      a <- xmlTreeParse(xml.file, useInternalNodes = TRUE)  
      
    ## Try to generate missing metadata:
    if(generate.missing == TRUE){
      message("Generating missing metadata...")
      ## Metadata template:
      if(md.type=="INSPIRE"){
        tmpl <- xmlTreeParse(system.file("INSPIRE_ISO19139.xml", package = "plotKML"), useInternalNodes = TRUE)
      } else {
        tmpl <- xmlTreeParse(system.file("FGDC.xml", package="plotKML"), useInternalNodes=TRUE)
      }
      top <- xmlRoot(tmpl)
      nx <- names(unlist(xmlToList(top)))
      ## compare the actual xml file and the template:
      cross <- compareXMLDocs(a=a, b=tmpl)
    
      ## Merge the existing Metadata file with existing template:
      if(length(cross[["inB"]])>0){
        for(i in 1:length(cross[["inB"]])){
          ## position of the missing node in the target doc:
          nodn <- attr(cross[["inB"]], "dimnames")[[1]][i]
          x_l <- strsplit(nx[grep(nodn, nx)], "\\.")[[1]] 
          ## TH: 7 levels - this is not the best implementation :(  -> it takes ca 3-5 seconds;
          for(j in 1:length(x_l)){
            if(j==1 & !x_l[j] %in% names(ml)) { ml <- append.XMLNode(ml, xmlNode(x_l[j], "")) }
            if(j==2 & !x_l[j] %in% names(ml[[x_l[1]]])) { ml[[x_l[1]]] <- append.XMLNode(ml[[x_l[1]]], xmlNode(x_l[j], "")) }
            if(j==3 & !x_l[j] %in% names(ml[[x_l[1]]][[x_l[2]]])) { ml[[x_l[1]]][[x_l[2]]] <- append.XMLNode(ml[[x_l[1]]][[x_l[2]]], xmlNode(x_l[j], "")) }
            if(j==4 & !x_l[j] %in% names(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]])) { ml[[x_l[1]]][[x_l[2]]][[x_l[3]]] <- append.XMLNode(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]], xmlNode(x_l[j], "")) }    
            if(j==5 & !x_l[j] %in% names(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]])) { ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]] <- append.XMLNode(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]], xmlNode(x_l[j], "")) } 
            if(j==6 & !x_l[j] %in% names(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]])) { ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]] <- append.XMLNode(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]], xmlNode(x_l[j], "")) }
            if(j==7 & !x_l[j] %in% names(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]][[x_l[6]]])) { ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]][[x_l[6]]] <- append.XMLNode(ml[[x_l[1]]][[x_l[2]]][[x_l[3]]][[x_l[4]]][[x_l[5]]][[x_l[6]]], xmlNode(x_l[j], "")) }
         }
      }}
      } 
       
      if(validate.schema == TRUE){
        if(md.type=="INSPIRE"){ 
           message("Connecting to the INSPIRE schema...")
           try(xsd <- xmlTreeParse(get("inspire_xsd", envir=plotKML.opts), isSchema=TRUE, useInternalNodes = TRUE))
        }
        if(md.type=="FGDC"){
           message("Connecting to the FGDC schema...")
           try(xsd <- xmlTreeParse(get("fgdc_xsd", envir=plotKML.opts), isSchema=TRUE, useInternalNodes = TRUE))        
        }
        ## validate if the schema is OK:
        try(xsd.val <- xmlSchemaValidate(xsd, a))
        ## print the results of validation:
        x <- NULL
        try(x <- sapply(xsd.val$errors, function(x){x["msg"][[1]]}))
        if(nzchar(x)){
          warning(paste(x), call. = FALSE)
        }
      }
      message(paste('Finished generating missing nodes. Open file \"', out.xml.file, '\" and add missing information.', sep=""))
    
    } else {
      warning(paste("Could not locate ", xml.file, ". See '?spMetadata' for more details.", sep=""), call.=FALSE)
    }
  }  ## If the metadata file does not exit, use the template available in the "/inst" directory:
  else {  
    ## Metadata template:
    if(md.type=="INSPIRE"){
      ret <- xmlTreeParse(system.file("INSPIRE_ISO19139.xml", package="plotKML"), useInternalNodes = TRUE)
    } 
    if(md.type=="FGDC"){
      ret <- xmlTreeParse(system.file("FGDC.xml", package="plotKML"), useInternalNodes = TRUE)
    }
    ml <- xmlRoot(ret) 
    
    ## Load static metadata entries (if they do not exist already):
    metadata.env(..., show.env = FALSE)
    ## Target variable: 
    Target_variable <- names(obj)[1] 
    ## Generate a name for the output XML file:
    if(missing(out.xml.file)){ 
      out.xml.file <- paste(normalizeFilename(deparse(substitute(obj, env=parent.frame()))), ".xml", sep="") 
    }
    ## Measurement resolution:
    Attribute_Measurement_Resolution <- get("Attribute_Measurement_Resolution", envir = metadata)
    if(Attribute_Measurement_Resolution==""){
      ## Estimate the numeric resolution by using the optimal number of bins in the histogram:
      if(is.numeric(obj@data[,Target_variable])){
         Attribute_Measurement_Resolution <- signif(diff(quantile(obj@data[,Target_variable], na.rm=TRUE, prob=c(.025,.975)))/(length(hist(obj@data[,Target_variable], breaks="FD", plot=FALSE)$breaks)*2), 2)
      }
    }
    
    ## The bounding box:
    if(md.type=="INSPIRE"){
      west <- get("Extent_West_Longitude", envir = metadata)
      east <- get("Extent_East_Longitude", envir = metadata)
      south <- get("Extent_South_Latitude", envir = metadata)
      north <- get("Extent_North_Latitude", envir = metadata)   
    } else {
      west <- get("West_Bounding_Coordinate", envir = metadata)
      east <- get("East_Bounding_Coordinate", envir = metadata)
      south <- get("South_Bounding_Coordinate", envir = metadata)
      north <- get("North_Bounding_Coordinate", envir = metadata)
    }
    if(any(list(west,east,south,north)=="")){
      message("Estimating the bounding box coordinates...")
      obj.ll <- reproject(obj)
      west <- min(coordinates(obj.ll)[,1])
      east <- max(coordinates(obj.ll)[,1])
      north <- max(coordinates(obj.ll)[,2])
      south <- min(coordinates(obj.ll)[,2])  
    }
    
    ## Automatically update some standard nodes (INSPIRE schema):  
    if(md.type=="INSPIRE"){
      ## Generate UUID:
      UUID <- get("UUID", envir = metadata)
      if(UUID==""){
        if(requireNamespace("uuid", quietly = TRUE)){
          xmlValue(ml[["fileIdentifier"]][[1]]) <- uuid::UUIDgenerate(use.time = FALSE)
        }
      }
      CI_RS_identifier <- get("CI_RS_identifier", envir = metadata)
      if(CI_RS_identifier==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["identifier"]][["RS_Identifier"]][["code"]][[1]]) <- uuid::UUIDgenerate(use.time = FALSE)
      }
      ## Responsible party:
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["organisationName"]][[1]]) <- get("MD_Organisation_name", envir = metadata)
      MD_Electronic_mail_address <- get("MD_Electronic_mail_address", envir = metadata) 
      if(!MD_Electronic_mail_address==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["contactInfo"]][["CI_Contact"]][["address"]][["CI_Address"]][["electronicMailAddress"]][[1]]) <- MD_Electronic_mail_address
      }
      MD_contact_URL <- get("MD_contact_URL", envir = metadata) 
      if(!MD_contact_URL==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["contactInfo"]][["CI_Contact"]][["onlineResource"]][["CI_OnlineResource"]][["linkage"]][[1]]) <- MD_contact_URL
      } 
      ## Citation title:
      CI_Citation_title <- get("CI_Citation_title", envir = metadata)
      if(CI_Citation_title==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["title"]][[1]]) <- normalizeFilename(deparse(substitute(obj)))
      } else {
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["title"]][[1]]) <- CI_Citation_title
      }
      CI_Unique_name <- get("CI_Unique_name", envir = metadata)
      if(!CI_Unique_name==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["identifier"]][["RS_Identifier"]][["code"]][[1]]) <- CI_Unique_name
      } 
      ## Abstract:
      MD_Abstract <- get("MD_Abstract", envir = metadata)
      if(!MD_Abstract==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["abstract"]][[1]]) <- MD_Abstract
      }
      ## Lineage:
      DQ_Lineage_statement <- get("DQ_Lineage_statement", envir = metadata)
      if(!DQ_Lineage_statement==""){
        xmlValue(ml[["dataQualityInfo"]][["DQ_DataQuality"]][["lineage"]][["LI_Lineage"]][["statement"]][[1]]) <- DQ_Lineage_statement
      }
      ## Use limitations:
      MD_Use_limitations <- get("MD_Use_limitations", envir = metadata)
      if(!MD_Use_limitations==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["resourceConstraints"]][["MD_Constraints"]][["useLimitation"]][[1]]) <- MD_Use_limitations
      }
      ## Scale number / resolution:
      MD_Equivalent_scale <- get("MD_Equivalent_scale", envir = metadata)
      xx <- which(names(ml[["identificationInfo"]][["MD_DataIdentification"]]) == "spatialResolution")  ## TH: there are two nodes of the same name!
      if(!MD_Equivalent_scale==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][[xx[1]]][["MD_Resolution"]][["equivalentScale"]][["MD_RepresentativeFraction"]][["denominator"]][["Integer"]]) <- MD_Equivalent_scale
      }
      MD_Resolution <- get("MD_Resolution", envir = metadata)
      if(!MD_Resolution==""){  
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][[xx[2]]][["MD_Resolution"]][["distance"]][["Distance"]]) <- MD_Resolution
        xmlAttrs(ml[["identificationInfo"]][["MD_DataIdentification"]][[xx[2]]][["MD_Resolution"]][["distance"]][["Distance"]])[[1]] <- "#m"
      }
      ## Coordinate system:
      MD_ReferenceSystem_Identifier <- get("MD_ReferenceSystem_Identifier", envir = metadata)
      if(MD_ReferenceSystem_Identifier==""){  
        xmlValue(ml[["referenceSystemInfo"]][["MD_ReferenceSystem"]][["referenceSystemIdentifier"]][["RS_Identifier"]][["code"]][[1]]) <- proj4string(obj)
        xmlValue(ml[["referenceSystemInfo"]][["MD_ReferenceSystem"]][["referenceSystemIdentifier"]][["RS_Identifier"]][["codeSpace"]][[1]]) <- "proj4"
      } else {
        xmlValue(ml[["referenceSystemInfo"]][["MD_ReferenceSystem"]][["referenceSystemIdentifier"]][["RS_Identifier"]][["code"]][[1]]) <- MD_ReferenceSystem_Identifier
        xmlValue(ml[["referenceSystemInfo"]][["MD_ReferenceSystem"]][["referenceSystemIdentifier"]][["RS_Identifier"]][["codeSpace"]][[1]]) <- get("MD_ReferenceSystem_type", envir = metadata)
      }
      
      ## Contact organization:
      xmlValue(ml[["contact"]][["CI_ResponsibleParty"]][["organisationName"]][[1]]) <- get("CI_Organisation_name", envir = metadata)
      xmlValue(ml[["contact"]][["CI_ResponsibleParty"]][["role"]][["CI_RoleCode"]][[1]]) <- get("CI_Role", envir = metadata)
      ## Author:
      CI_Electronic_mail_address <- get("CI_Electronic_mail_address", envir = metadata)
      if(!CI_Electronic_mail_address==""){
        xmlValue(ml[["contact"]][["CI_ResponsibleParty"]][["contactInfo"]][["CI_Contact"]][["address"]][["CI_Address"]][["electronicMailAddress"]][[1]]) <- CI_Electronic_mail_address
      }
      
      ## Thesaurus:
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["descriptiveKeywords"]][["MD_Keywords"]][["thesaurusName"]][["CI_Citation"]][["title"]][[1]]) <- get("MD_Thesaurus_name", envir = metadata)
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["descriptiveKeywords"]][["MD_Keywords"]][["thesaurusName"]][["CI_Citation"]][["date"]][["CI_Date"]][["date"]][["Date"]]) <- get("MD_Thesaurus_date", envir = metadata)
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["descriptiveKeywords"]][["MD_Keywords"]][["thesaurusName"]][["CI_Citation"]][["date"]][["CI_Date"]][["dateType"]][["CI_DateTypeCode"]]) <- get("MD_Thesaurus_date_type", envir = metadata)
      MD_Keyword <- get("MD_Keyword", envir = metadata)
      if(!MD_Keyword==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["descriptiveKeywords"]][["MD_Keywords"]][["keyword"]][[1]]) <- MD_Keyword
      }
      
      ## Language:
      xmlValue(ml[["language"]][["LanguageCode"]]) <- get("Language_code", envir = metadata)
      xmlAttrs(ml[["language"]][["LanguageCode"]])[[2]] <- get("Language_code", envir = metadata)
      xmlValue(ml[["characterSet"]][["MD_CharacterSetCode"]][[1]]) <- get("MD_Character_set_code", envir = metadata)
      xmlValue(ml[["hierarchyLevel"]][["MD_ScopeCode"]][[1]]) <- get("MD_scope_code", envir = metadata)
      ## Citation date:
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["citation"]][["CI_Citation"]][["date"]][["CI_Date"]][["date"]][["Date"]]) <- get("CI_Citation_date", envir = metadata)
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["language"]][["LanguageCode"]]) <- get("MD_Language_code", envir = metadata)
      xmlAttrs(ml[["identificationInfo"]][["MD_DataIdentification"]][["language"]][["LanguageCode"]])[[2]] <- get("MD_Language_code", envir = metadata)
      ## All topic categories:
      MD_Topic_category_code <- get("MD_Topic_category_code", envir = metadata)
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["topicCategory"]][["MD_TopicCategoryCode"]]) <- MD_Topic_category_code[1]
      if(length(MD_Topic_category_code) > 1) {
        xi <- length(MD_Topic_category_code) - 1
        xt <- as.vector(which(names(ml[["identificationInfo"]][["MD_DataIdentification"]]) == "topicCategory")[1])
        xm <- xmlSize(ml[["identificationInfo"]][["MD_DataIdentification"]])
        j <- 0
        for (i in (xt + 1):xm) {
          ml[["identificationInfo"]][["MD_DataIdentification"]][xm + xi + j][[1]] <- ml[["identificationInfo"]][["MD_DataIdentification"]][i][[1]]
          j <- j + 1
        }
        for(i in 2:length(MD_Topic_category_code)) {
          ml[["identificationInfo"]][["MD_DataIdentification"]][xt + i - 1][[1]] <- ml[["identificationInfo"]][["MD_DataIdentification"]][["topicCategory"]]
          xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][xt + i - 1][[1]][["MD_TopicCategoryCode"]]) <- MD_Topic_category_code[i]
       }
      }
      ## Data quality info:
      DQ_Citation_title <- get("DQ_Citation_title", envir = metadata)
      if(!DQ_Citation_title==""){
        xmlValue(ml[["dataQualityInfo"]][["DQ_DataQuality"]][["report"]][["DQ_DomainConsistency"]][["result"]][["DQ_ConformanceResult"]][["specification"]][["CI_Citation"]][["title"]][[1]]) <- DQ_Citation_title
      }
      DQ_Citation_date <- get("DQ_Citation_date", envir = metadata)
      if(!DQ_Citation_title==""){
        xmlValue(ml[["dataQualityInfo"]][["DQ_DataQuality"]][["report"]][["DQ_DomainConsistency"]][["result"]][["DQ_ConformanceResult"]][["specification"]][["CI_Citation"]][["date"]][["CI_Date"]][["date"]][["Date"]]) <- DQ_Citation_date
      }
      xmlValue(ml[["dataQualityInfo"]][["DQ_DataQuality"]][["report"]][["DQ_DomainConsistency"]][["result"]][["DQ_ConformanceResult"]][["specification"]][["CI_Citation"]][["date"]][["CI_Date"]][["dateType"]][["CI_DateTypeCode"]]) <- get("DQ_Citation_date_type", envir = metadata)
     ## Resource locator...
     CI_Online_resource_URL <- get("CI_Online_resource_URL", envir = metadata) 
     if(!CI_Online_resource_URL==""){
        xmlValue(ml[["distributionInfo"]][["MD_Distribution"]][["transferOptions"]][["MD_DigitalTransferOptions"]][["onLine"]][["CI_OnlineResource"]][["linkage"]][["URL"]]) <- CI_Online_resource_URL
     }
     if(validate.schema==TRUE){
       try.connection <- try(url(CI_Online_resource_URL, open = 'rb'))
       try.error <- inherits(try.connection, "try-error")
       if(try.error == T){
         warning("'CI_Online_resource_URL' points to an invalid URL.")
       }
     }
     ## Bounding box:
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["geographicElement"]][["EX_GeographicBoundingBox"]][["westBoundLongitude"]][["Decimal"]]) <- west
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["geographicElement"]][["EX_GeographicBoundingBox"]][["eastBoundLongitude"]][["Decimal"]]) <- east
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["geographicElement"]][["EX_GeographicBoundingBox"]][["northBoundLatitude"]][["Decimal"]]) <- north
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["geographicElement"]][["EX_GeographicBoundingBox"]][["southBoundLatitude"]][["Decimal"]]) <- south
      ## Dates:
      Time_period_begin <- get("Time_period_begin", envir = metadata) 
      if(!Time_period_begin==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["temporalElement"]][["EX_TemporalExtent"]][["extent"]][["TimePeriod"]][["beginPosition"]]) <- as(Time_period_begin, "character")
      }
      Time_period_end <- get("Time_period_end", envir = metadata) 
      if(!Time_period_end==""){
        xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["extent"]][["EX_Extent"]][["temporalElement"]][["EX_TemporalExtent"]][["extent"]][["TimePeriod"]][["endPosition"]]) <- as(Time_period_end, "character")
      }
      ## Metadata standard:
      xmlValue(ml[["metadataStandardName"]][[1]]) <- get("MD_Standard_name", envir = metadata)
      xmlValue(ml[["metadataStandardVersion"]][[1]]) <- get("MD_Standard_version", envir = metadata)
      ## Contact:
      xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["role"]][["CI_RoleCode"]]) <- get("MD_Role", envir = metadata)
     xmlAttrs(ml[["identificationInfo"]][["MD_DataIdentification"]][["pointOfContact"]][["CI_ResponsibleParty"]][["role"]][["CI_RoleCode"]])[[2]] <- get("MD_Role", envir = metadata)
     ## Use constraints:
     xx <- which(names(ml[["identificationInfo"]][["MD_DataIdentification"]]) == "resourceConstraints")  ## TH: there are two nodes of the same name!
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][[xx[2]]][["MD_LegalConstraints"]][["accessConstraints"]][["MD_RestrictionCode"]]) <- get("MD_Access_constraints", envir = metadata)
     xmlValue(ml[["identificationInfo"]][["MD_DataIdentification"]][[xx[2]]][["MD_LegalConstraints"]][["otherConstraints"]][[1]]) <- get("MD_Other_restrictions", envir = metadata)
  
    }
             
    ## Automatically update some standard nodes (FGDC schema):   
    if(md.type=="FGDC"){
      ## Location name:
      Indirect_Spatial_Reference <- get("Indirect_Spatial_Reference", envir = metadata)
      if(Indirect_Spatial_Reference=="" & GoogleGeocode == TRUE){
        if(requireNamespace("rjson", quietly = TRUE)){
          googleurl <- url(paste("http://maps.googleapis.com/maps/api/geocode/json?latlng=",  round(mean(obj.ll@bbox[2,]),3), ",", round(mean(obj.ll@bbox[1,]),3), "&sensor=false", sep=""))
          try(Indirect_Spatial_Reference <- rjson::fromJSON(file=googleurl)[["results"]][[1]][["formatted_address"]])
          close(googleurl)
        }
      }

      ## Citation title:
      Citation_title <- get("Citation_title", envir = metadata)
      if(Citation_title==""){
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["title"]]) <- normalizeFilename(deparse(substitute(obj)))
      } else {
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["title"]]) <- Citation_title
      }
      ## Data format:
      xmlValue(ml[["distinfo"]][["stdorder"]][["digform"]][["digtinfo"]][["formcont"]]) <- class(obj)[1]
      ## Data citation:
      xmlValue(ml[["eainfo"]][["overview"]][["eadetcit"]]) <- paste("http://CRAN.R-project.org/package=", attr(class(obj), "package"), "/", sep="")
      ## Entity type:
      xmlValue(ml[["eainfo"]][["detailed"]][["enttyp"]][["enttypl"]]) <- class(obj@data[, Target_variable])
      ## Range of values:
      if(is.numeric(obj@data[,Target_variable])){
        xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdomv"]][["rdom"]][["rdommin"]]) <- min(obj@data[,Target_variable], na.rm=TRUE) 
        xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdomv"]][["rdom"]][["rdommax"]]) <- max(obj@data[,Target_variable], na.rm=TRUE)
      }  
      ## Attribute label:
      xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrlabl"]]) <- Target_variable
      Attribute_Definition <- get("Attribute_Definition", envir = metadata)
      if(!Attribute_Definition==""){
        xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdef"]]) <- Attribute_Definition
      }
      ## MB: geoform node:
      Geospatial_Data_Presentation_Form <- get("Geospatial_Data_Presentation_Form", envir = metadata)
      if(!Geospatial_Data_Presentation_Form==""){
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["geoform"]]) <- Geospatial_Data_Presentation_Form
      }
      ## MB: 'pubinfo' nodes:
      Citation_Publisher <- get("Citation_Publisher", envir = metadata)
      if(!Citation_Publisher==""){
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["pubinfo"]][["publish"]]) <- Citation_Publisher
      }
      Publication_Place <- get("Publication_Place", envir = metadata)
      if(!Publication_Place==""){
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["pubinfo"]][["pubplace"]]) <- Publication_Place
      }
      Citation_URL <- get("Citation_URL", envir = metadata)
      if(!Citation_URL==""){
        xmlValue(ml[["idinfo"]][["citation"]][["citeinfo"]][["onlink"]]) <- Citation_URL
      }
      ## MB: insert 'edition' node:
      Purpose <- get("Purpose", envir = metadata)
      if(!Purpose==""){
        xmlValue(ml[["idinfo"]][["descript"]][["purpose"]]) <- Purpose
      }
      Abstract <- get("Abstract", envir = metadata)
      if(!Abstract==""){
        xmlValue(ml[["idinfo"]][["descript"]][["abstract"]]) <- Abstract
      }
      Supplemental_Information <- get("Supplemental_Information", envir = metadata)    
      if(!Supplemental_Information==""){
        xmlValue(ml[['idinfo']][['descript']]['supplinf']) <- Supplemental_Information
      }
      ## temporal extent:
      Beginning_Date <- get("Beginning_Date", envir = metadata)
      if(!Beginning_Date==""){
        xmlValue(ml[["idinfo"]][["timeperd"]][["timeinfo"]][["rngdates"]][["begdate"]]) <- Beginning_Date
      }
      Ending_Date <- get("Ending_Date", envir = metadata)
      if(!Ending_Date==""){
        xmlValue(ml[["idinfo"]][["timeperd"]][["timeinfo"]][["rngdates"]][["enddate"]]) <- Ending_Date
      }
      ## keywords theme:
      Theme_Keyword <- get("Theme_Keyword", envir = metadata)
      if(!Theme_Keyword==""){
        for(i in 1:length(Theme_Keyword)) {
          if(i != length(Theme_Keyword)) {
            ml[['idinfo']][['keywords']][['theme']][i+2] <- ml[['idinfo']][['keywords']][['theme']][i+1]
          }
          xmlValue(ml[['idinfo']][['keywords']][['theme']][i+1][['themekey']]) <- Theme_Keyword[i]
       }
      }
      Theme_Keyword_Thesaurus <- get("Theme_Keyword_Thesaurus", envir = metadata)    
      if(!Theme_Keyword_Thesaurus==""){
        for(i in 1:length(Theme_Keyword_Thesaurus)){
          if(i != length(Theme_Keyword_Thesaurus)) {
            ml[['idinfo']][['keywords']][['place']][i+2] <- ml[['idinfo']][['keywords']][['place']][i+1]
          } 
          xmlValue(ml[['idinfo']][['keywords']][['place']][i+1][['placekey']]) <- Theme_Keyword_Thesaurus[i]
        }
      } 
      ## Measurement resolution:
      xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdomv"]][["rdom"]][["attrmres"]]) <- Attribute_Measurement_Resolution    
      Attribute_Units_of_Measure <- get("Attribute_Units_of_Measure", envir = metadata)
      if(!Attribute_Units_of_Measure==""){
        xmlValue(ml[["eainfo"]][["detailed"]][["attr"]][["attrdomv"]][["rdom"]][["attrunit"]]) <- Attribute_Units_of_Measure
      }
      ## Estimate object size:
      xmlValue(ml[["distinfo"]][["stdorder"]][["digform"]][["digtinfo"]][["transize"]]) <-  paste(signif(object.size(obj@data[, Target_variable])/1024, 4), "Kb (uncompressed)") 
      ## MB: insert status information (progress, update):
      Status_Progress <- get("Status_Progress", envir = metadata)
      if(length(Status_Progress)==1L){
        xmlValue(ml[["idinfo"]][["status"]][["progress"]]) <- Status_Progress
      }
      Maintenance_and_Update_Frequency <- get("Maintenance_and_Update_Frequency", envir = metadata)  
      if(!Maintenance_and_Update_Frequency==""){
        xmlValue(ml[["idinfo"]][["status"]][["update"]]) <- Maintenance_and_Update_Frequency
      } 
      ## Bounding box:
      xmlValue(ml[["idinfo"]][["spdom"]][["bounding"]][["westbc"]]) <- west
      xmlValue(ml[["idinfo"]][["spdom"]][["bounding"]][["eastbc"]]) <- east
      xmlValue(ml[["idinfo"]][["spdom"]][["bounding"]][["northbc"]]) <- north
      xmlValue(ml[["idinfo"]][["spdom"]][["bounding"]][["southbc"]]) <- south  
      xmlValue(ml[["spref"]][["horizsys"]][["geograph"]][["geogunit"]]) <- "Decimal degrees"
      ## Representation type:
      Direct_Spatial_Reference_Method <- get("Direct_Spatial_Reference_Method", envir = metadata)[1]
      if(!Direct_Spatial_Reference_Method==""){
        xmlValue(ml[['spdoinfo']][['direct']]) <- Direct_Spatial_Reference_Method
      }
      ## metadata generation time:    
      xmlValue(ml[["metainfo"]][["metd"]]) <- Sys.Date()
      ## Spatial reference:
      if(!Indirect_Spatial_Reference==""){
        xmlValue(ml[["spdoinfo"]][["indspref"]]) <- Indirect_Spatial_Reference
      }
      ## number of elements:
      xmlValue(ml[["spdoinfo"]][["ptvctinf"]][["sdtsterm"]][["ptvctcnt"]]) <- nrow(obj@data)
      ## System information:
      xmlValue(ml[["idinfo"]][["native"]]) <- get("Native_Data_Set_Environment", envir = metadata)
      xmlValue(ml[["idinfo"]][["ptcontac"]][["cntinfo"]][["cntorgp"]][["cntper"]]) <- get("Contact_Information_Person", envir = metadata)
      Metadata_Contact_Position <- get("Metadata_Contact_Position", envir = metadata)
      if(!Metadata_Contact_Position==""){
        xmlValue(ml[["metainfo"]][["metc"]][["cntinfo"]][["cntpos"]]) <- get("Metadata_Contact_Position", envir = metadata)
      }
    }

  }
  message("Finished generating standard metadata...\nManual editing of the metadata file might be required.")
  
  ## attach friendly metadata names:
  ny <- unlist(xmlToList(ml))
  ## convert to a table:
  met <- data.frame(metadata=names(ny), value=paste(ny))
  ## add friendly names:
  mdnames <- read.csv(system.file("mdnames.csv", package="plotKML"))
  field_names <- merge(met, mdnames[,c("metadata","field.names")], by="metadata", all.x=TRUE, all.y=FALSE, sort=FALSE)[,"field.names"]
  field_names <- field_names[!is.na(field_names)]
  ## generate metadata doc:
  saveXML(ml, out.xml.file)
  doc = xmlInternalTreeParse(out.xml.file)
  
  ## color palette:
  if(missing(bounds)){
    if(is.numeric(obj@data[,Target_variable])){
      bounds <- seq(range(obj@data[,Target_variable], na.rm = TRUE, finite = TRUE)[1], range(obj@data[,Target_variable], na.rm = TRUE, finite = TRUE)[2], Attribute_Measurement_Resolution/2)  ## half the numeric resolution!
      bounds.c <- signif((bounds[-1]+bounds[-length(bounds)])/2, signif.digit)
      if(missing(legend_names)) { legend_names <- as.character(bounds.c) }
    } 
    else { 
      x <- as.factor(obj@data[,Target_variable])
      bounds <- c(0, seq_along(levels(x)))
      bounds.c <- bounds[-1]
      if(missing(legend_names)) { legend_names <- as.character(levels(x)) }   
    } }
  
  ## generate a palette:
  if(missing(colour_scale)){
    if(is.numeric(obj@data[,Target_variable])){
      colour_scale <- get("colour_scale_numeric", envir = plotKML.opts)
    }
    else { 
      colour_scale <- get("colour_scale_factor", envir = plotKML.opts) 
    }
  }
  
  if(missing(icons)){ 
    icons <- rep("", length(legend_names))
  }
  
  cols <- colorRamp(colour_scale, space="rgb", interpolate = "linear")
  cdata <- scales::rescale(bounds.c)
  if(missing(color)){
    color <- rgb(cols(cdata)/255)
  }
  
  ## make a spatial palette:
  pal <- new("sp.palette", type=class(obj@data[,Target_variable]), bounds=bounds, color = color, names = legend_names, icons = icons)
  ## make a SpatialMetadata object:
  spmd <- new("SpatialMetadata",  xml=doc, field.names=paste(field_names), palette=pal, sp=as(obj, "Spatial")) 
  return(spmd)
  
}

setMethod("spMetadata", "Spatial", .spMetadata.Spatial)

## Raster objects:
.spMetadata.Raster <- function(obj, bounds = NULL, color = NULL, ...){
  
  if(!is.null(bounds)) {bounds <- obj@legend@values}
  if(!is.null(color)) {color <- obj@legend@color}
  
  ## convert a Raster layer to SGDF:
  if(nlayers(obj) > 1){
    obj <- raster(obj, layer = 1)
  }
  obj <- as(obj, "SpatialGridDataFrame") 
  
  .spMetadata.Spatial(obj, ...)
}

setMethod("spMetadata", "RasterLayer", .spMetadata.Raster)


################## READ METADATA ##############

## Read metadata from a xml.file and convert to a table:
read.metadata <- function(xml.file, delim.sign, full.names){
  
  if(missing(full.names)){    
    full.names = read.csv(system.file("mdnames.csv", package="plotKML"))      
  }
  ret <- xmlTreeParse(xml.file, useInternalNodes = TRUE)
  top <- xmlRoot(ret)
  nx <- unlist(xmlToList(top))  ## , addAttributes=FALSE
  ## convert to a table:
  if(!missing(delim.sign)){
    nxn <- gsub("\\.", delim.sign, attr(nx, "names"))
  } else {
    nxn <- attr(nx, "names")
  }
  met <- data.frame(metadata = nxn, value=paste(nx), stringsAsFactors = FALSE)
  ## add more friendly names:
  metm <- merge(x=met, y=full.names[,c("metadata","field.names")], by="metadata", all.x=TRUE, sort=FALSE)
  
  return(metm[,c(1,3,2)])
}

# end of script;