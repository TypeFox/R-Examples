CDA2ODM <- function( CDAfile="" ) 
{

if (CDAfile == "") CDAfile <- file.choose()
currentdir <- getwd()
setwd(dirname(CDAfile))

CDA = xmlRoot(xmlTreeParse(CDAfile))

filename <- paste(CDAfile, ".odm.xml", sep="")
print(paste("Generating", filename))
sink(filename)

Description <- paste( xmlValue(CDA[["title"]]), xmlAttrs(CDA[["templateId"]])["assigningAuthorityName"], sep=" ")

cat("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\r", 
    "<ODM xmlns=\"http://www.cdisc.org/ns/odm/v1.3\"\r", 
    "Description=\"", Description, "\"\r",
    "ODMVersion=\"1.3.1\" ",
    "CreationDateTime=\"", format(Sys.time(), "%Y-%m-%dT%X"), "+01:00", "\"\r",
    "FileOID=\"", filename, "\"\r",
    "FileType=\"Snapshot\"", ">\r\r",
sep="")

StudyOID <- paste( xmlAttrs(CDA[["id"]])["root"], ".", xmlAttrs(CDA[["id"]])["extension"], sep="" )
StudyName <- paste( xmlValue(CDA[["title"]]), xmlAttrs(CDA[["templateId"]])["assigningAuthorityName"], sep=" ")
StudyDescription <-  paste( xmlAttrs(CDA[["code"]])["displayName"], xmlAttrs(CDA[["code"]])["codeSystemName"], xmlAttrs(CDA[["code"]])["code"],  
                            xmlAttrs(CDA[["id"]])["assigningAuthorityName"], sep=" ")
ProtocolName <- paste( xmlAttrs(CDA[["templateId"]])["root"], xmlAttrs(CDA[["templateId"]])["assigningAuthorityName"], sep=" ")
cat("<Study OID=\"", StudyOID, "\">\r", 
    "<GlobalVariables>\r",
    "<StudyName>", StudyName, "</StudyName>\r",
    "<StudyDescription>", StudyDescription, "</StudyDescription>\r",
    "<ProtocolName>", ProtocolName, "</ProtocolName>\r",
    "</GlobalVariables>\r", 
sep="")

cat("<BasicDefinitions>\r<MeasurementUnit OID=\"MU.1\" Name=\"unit\">\r",
    "<Symbol><TranslatedText xml:lang=\"en\">Unit</TranslatedText></Symbol>\r",
    "</MeasurementUnit>\r</BasicDefinitions>\r",
sep="")

cat("<MetaDataVersion OID=\"MD.1\" Name=\"Metadataversion\">\r",
    "<Protocol><StudyEventRef StudyEventOID=\"SE.1\" OrderNumber=\"1\" Mandatory=\"Yes\" />", 
    "</Protocol>\r",
    "<StudyEventDef OID=\"SE.1\" Name=\"", xmlValue(CDA[["title"]]), "\" Repeating=\"No\" Type=\"Unscheduled\">\r",
    "<FormRef FormOID=\"F.1\" OrderNumber=\"1\" Mandatory=\"No\" />\r",
    "</StudyEventDef>\r",
sep="")

cat("<FormDef OID=\"F.1\" Name=\"", xmlValue(CDA[["title"]]), "\" Repeating=\"No\">\r",
    "   <ItemGroupRef ItemGroupOID=\"IG.1\" Mandatory=\"Yes\" OrderNumber=\"1\" />\r", 
    "</FormDef>\r", 
    "<ItemGroupDef OID=\"IG.1\" Name=\"ClinicalDocument\" Repeating=\"No\">\r",
sep="")

# extract and normalize itemnames
itemnames <- NULL
for (i in 1:length(CDA))
{
   itemnames <- c( itemnames, names(unlist(CDA[i])) )
}
itemnames <- gsub(".children","",itemnames, fixed=T)
itemnames <- gsub(".text.value$","",itemnames, perl=T)
itemnames <- itemnames[!grepl(".namespace$",itemnames, perl=T)]
itemnames <- itemnames[substr(itemnames, nchar(itemnames) - 4,nchar(itemnames)) != ".name"]
itemnames <- itemnames[substr(itemnames, nchar(itemnames) - 12,nchar(itemnames)) != "comment.value"]

imax <- length(itemnames)

for (i in 1:imax)
{
   cat("   <ItemRef ItemOID=\"I.", i, "\" Mandatory=\"Yes\" OrderNumber=\"", i, "\" />\r", sep="")
}
cat("</ItemGroupDef>\r")

for (i in 1:imax)
{
   tmp <- unlist(strsplit(itemnames[i],".", fixed=T))
   last <- length(tmp)
   if (last < 3) 
   { 
      shortname <- tmp[last] 
   }
   else
   {
      if (tmp[last -1] == "attributes") { shortname <- paste(tmp[last -2], tmp[last], sep=".") }
      else { shortname <- tmp[last] }
   }
   cat("<ItemDef OID=\"I.", i, "\" Name=\"", itemnames[i], 
       "\" DataType=\"string\">\r", 
       "   <Question>\r",
       "      <TranslatedText xml:lang=\"en\">", shortname, "</TranslatedText>\r",
       "   </Question>\r", 
       "</ItemDef>\r", sep="")
}

cat("</MetaDataVersion>\r")
cat("</Study>\r\r")

cat("<AdminData>\r")
cat("   <User OID=\"USR.1\">\r")
cat("      <LoginName>CDA2ODM</LoginName>\r", sep="") 
cat("      <DisplayName />\r") 
cat("      <FullName>CDA2ODM converter</FullName>\r", sep="") 
cat("      <FirstName>NN</FirstName>\r", sep="")
cat("      <LastName>NN</LastName>\r", sep="") 
cat("      <Organization>NN</Organization>\r", sep="")
cat("      <Pager />\r")
cat("   </User>\r")
cat("   <Location OID=\"LOC.1\" Name=\"NN\" LocationType=\"Other\">\r", sep="")
cat("      <MetaDataVersionRef StudyOID=\"", StudyOID, "\" MetaDataVersionOID=\"MD.1\" ",
           "EffectiveDate=\"", format(Sys.time(), "%Y-%m-%d"), "\" />\r", sep="") 
cat("   </Location>\r")
cat("</AdminData>\r")
cat("</ODM>\r")

sink(); print("Finished\r")
}
