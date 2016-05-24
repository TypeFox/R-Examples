# Martin Dugas 2013
# license: GPL

R2ODM <- function(odmdata) {

StudyOID  <- attr(odmdata,"StudyOID")
if (is.null(StudyOID)) StudyOID <- "StudyOID"
Sponsor   <- attr(odmdata,"Sponsor")
if (is.null(Sponsor)) Sponsor <- ""
Condition <- attr(odmdata,"Condition")
if (is.null(Condition)) Condition <- ""
StudyName <- attr(odmdata,"StudyName")
if (is.null(StudyName)) StudyName <- "StudyName"
StudyDescription <- attr(odmdata,"StudyDescription")
if (is.null(StudyDescription)) StudyDescription <- ""
Form <- attr(odmdata,"Form")
if (is.null(Form)) Form <- "Form1"
FirstName <- attr(odmdata,"FirstName")
if (is.null(FirstName)) FirstName <- ""
LastName <-  attr(odmdata,"LastName")
if (is.null(LastName)) LastName <- ""
Organization <- attr(odmdata,"Organization")
if (is.null(Organization)) Organization <- "DefaultOrganization"

cat("R2ODM\n")
filename <- Form
if (StudyOID != "") filename <- paste(filename, StudyOID, sep="_")
filename <- paste(filename, ".xml", sep="")
cat(paste("Generating \"", filename, "\"\n",sep=""))

sink(filename)
cat("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\r", 
    "<ODM xmlns=\"http://www.cdisc.org/ns/odm/v1.3\"\r", 
    "Description=\"Sponsor=", Sponsor, ",Condition=", Condition, ",Form=", Form, "\"\r",
    "ODMVersion=\"1.3.1\" ",
    "CreationDateTime=\"", format(Sys.time(), "%Y-%m-%dT%X"), "+01:00", "\"\r",
    "FileOID=\"", filename, "\"\r",
    "FileType=\"Snapshot\"", ">\r",
sep="")

cat("<Study OID=\"", StudyOID, "\">\r", 
    "<GlobalVariables>\r",
    "<StudyName>", StudyName, "</StudyName>\r",
    "<StudyDescription>", StudyDescription, "</StudyDescription>\r",
    "<ProtocolName>", Sponsor, "</ProtocolName>\r",
    "</GlobalVariables>\r", 
sep="")

cat("<BasicDefinitions>\r<MeasurementUnit OID=\"MU.1\" Name=\"unit\">\r",
    "<Symbol><TranslatedText xml:lang=\"en\">Unit</TranslatedText></Symbol>\r",
    "</MeasurementUnit>\r</BasicDefinitions>\r",
sep="")

cat("<MetaDataVersion OID=\"MD.1\" Name=\"Metadataversion\">\r",
    "<Protocol><StudyEventRef StudyEventOID=\"SE.1\" OrderNumber=\"1\" Mandatory=\"Yes\" />", 
    "</Protocol>\r",
    "<StudyEventDef OID=\"SE.1\" Name=\"", Condition, "\" Repeating=\"No\" Type=\"Unscheduled\">\r",
    "<FormRef FormOID=\"F.1\" OrderNumber=\"1\" Mandatory=\"No\" />\r",
    "</StudyEventDef>\r",
sep="")

# FormDef with ItemGroupRef
n <- length(names(odmdata))
cat("<FormDef OID=\"F.1\" Name=\"", Form, "\" Repeating=\"No\">\r",sep="")
IGU <- attr(odmdata,"itemgroups_unique")
if (is.null(IGU)) IGU <- as.vector("IG.1")
for (i in 1:length(IGU)) {
   cat("<ItemGroupRef ItemGroupOID=\"", IGU[i], "\" Mandatory=\"Yes\" OrderNumber=\"",i,"\" />\r", sep="")
}
cat("</FormDef>\r")

## itemgroup und itemdefs
IG <- attr(odmdata,"itemgroups")
if (is.null(IG)) IG <- rep("IG.1",n)
IGNU <- attr(odmdata,"itemgroupnames_unique")
IGNU2 <- attr(odmdata,"itemgroupnames_unique_en")
IGNU3 <- attr(odmdata,"itemgroupnames_unique_de")
IGAlias <- attr(odmdata,"IGAlias")
mx <- attr(odmdata,"Alias_Items")
Contexts <- unique(mx[,2])

if (is.null(IGNU)) IGNU <- as.vector("IG.1")
for (i in 1:length(IGU)) {
   itemno <- 0
   cat("<ItemGroupDef OID=\"", IGU[i], "\" Name=\"", IGNU[i], "\" Repeating=\"No\">\n",sep="")
   if (!is.null(IGNU2))
   {
      cat("   <Description>", "\n", "      <TranslatedText xml:lang=\"en\">", IGNU2[i], "</TranslatedText>\n",sep="")
      if (!is.null(IGNU3))
      {
         cat("      <TranslatedText xml:lang=\"de\">", IGNU3[i], "</TranslatedText>\n",sep="")
      }
	  cat("   </Description>\n")
   }
   for (ii in 1:n) {
      if ( IGU[i] == IG[ii]) {
         itemno <- itemno + 1
         cat("   <ItemRef ItemOID=\"", names(odmdata)[ii], 
             "\" Mandatory=\"Yes\" OrderNumber=\"", itemno, "\" />\n", sep="")
      }
   }
   if (length(IGAlias) > 0) for (ii in seq(1,length(IGAlias),3)) 
   {
      if (IGAlias[ii] == IGU[i])
	  {
	     cat("   <Alias Context=\"", IGAlias[ii+1], "\" Name=\"", IGAlias[ii+2], "\" />\n", sep="")
	  }
   }
   
   cat("</ItemGroupDef>\n")
}
# ItemDefs
varlabels <- attr(odmdata,"varlabels")
if (is.null(varlabels)) varlabels <- names(odmdata)
varlabels_en <- attr(odmdata,"varlabels_en")
if (is.null(varlabels_en)) varlabels_en <- varlabels
varlabels_de <- attr(odmdata,"varlabels_de")
if (is.null(varlabels_de)) varlabels_de <- varlabels
f3 <- function (x) gsub("<","&lt;",x)
f4 <- function (x) gsub(">","&gt;",x)
varlabels <- sapply(varlabels,f3,USE.NAMES=F)
varlabels <- sapply(varlabels,f4,USE.NAMES=F)
varlabels_en <- sapply(varlabels_en,f3,USE.NAMES=F)
varlabels_en <- sapply(varlabels_en,f4,USE.NAMES=F)
varlabels_de <- sapply(varlabels_de,f3,USE.NAMES=F)
varlabels_de <- sapply(varlabels_de,f4,USE.NAMES=F)

DataTypes <- attr(odmdata,"DataTypes")
if (is.null(DataTypes)) {
   DataTypes <- rep("string",n)
   for (i in 1:n) {
      if (typeof(odmdata[,i]) == "integer") DataTypes[i] <- "integer"
      if (typeof(odmdata[,i]) == "double")  DataTypes[i] <- "float"
   }
}
CL <- 0
CLxml <- ""
for (i in 1:n)
   {
   cat("<ItemDef OID=\"", names(odmdata)[i], "\" Name=\"", varlabels[i], "\" DataType=\"",DataTypes[i], "\">\n", sep="")
   cat("   <Question>\n",
       "      <TranslatedText xml:lang=\"en\">", varlabels_en[i], "</TranslatedText>\n",
       "      <TranslatedText xml:lang=\"de\">", varlabels_de[i], "</TranslatedText>\n",
       "   </Question>\n", sep="")
   # codelist
   if (nlevels(odmdata[,i]) >0) {
      # codelist
      itemcontexts <- attr(odmdata[,i],"Alias_Contexts")
      CL <- CL + 1
      cat("   <CodeListRef CodeListOID=\"CL.", CL, "\"/>\n", sep="")
      CLxml <- paste(CLxml, 
               "<CodeList OID=\"CL.", CL, "\" Name=\"", varlabels[i], "\" DataType=\"string\">\n",
                  sep="")
      if (is.null(attr(odmdata[,i],"labels"))) attr(odmdata[,i],"labels") <- levels(odmdata[,i]) 
      if (is.null(attr(odmdata[,i],"labels2"))) attr(odmdata[,i],"labels2") <- attr(odmdata[,i],"labels") 
      for (ii in 1: nlevels(odmdata[,i]) ) {
         CLxml <- paste(CLxml,
                  "   <CodeListItem CodedValue=\"", levels(odmdata[,i])[ii], "\">\n",
                  "      <Decode>\n",
                  "         <TranslatedText xml:lang=\"en\">", attr(odmdata[,i],"labels")[ii], "</TranslatedText>\n",
                  "         <TranslatedText xml:lang=\"de\">", attr(odmdata[,i],"labels2")[ii], "</TranslatedText>\n",
                  "      </Decode>\n", sep="")
         if (length(itemcontexts > 0)) for (iii in 1:length(itemcontexts))
         {
            aname <- attr(odmdata[,i],itemcontexts[iii])[ii]
            if (aname != "") {
               CLxml <- paste(CLxml, "      <Alias Context=\"", itemcontexts[iii], "\" Name=\"", aname, "\" />\n", sep="")
            }
         }
         CLxml <- paste(CLxml, "   </CodeListItem>\n", sep="")
      }
      CLxml <- paste(CLxml, "</CodeList>\n", sep="")
   }

   # Alias for items
   for (ii in 1:nrow(mx))
   {
      if (mx[ii,1] == names(odmdata)[i])
      {
         cat("   <Alias Context=\"", mx[ii,2], "\" Name=\"", mx[ii,3], "\" />\n", sep="")
      }
   }

   #if (length(Contexts) > 0) for (ii in 1:length(Contexts))
   #{
   #   aname <- attr(odmdata,Contexts[ii])[i]
   #   if (aname == "") next
   #   cat("   <Alias Context=\"", Contexts[ii], 
   #          "\" Name=\"", aname, "\" />\n", sep="")
   #}

   cat("</ItemDef>\n")
}
cat(CLxml)

cat("</MetaDataVersion>\n")
cat("</Study>\n\n")
cat("<AdminData>\n")
cat("   <User OID=\"USR.1\">\n")
cat("      <LoginName>", LastName, "</LoginName>\n", sep="") 
cat("      <DisplayName />\n") 
cat("      <FullName>", FirstName, " ", LastName, "</FullName>\n", sep="") 
cat("      <FirstName>", FirstName, "</FirstName>\n", sep="")
cat("      <LastName>", LastName, "</LastName>\n", sep="") 
cat("      <Organization>", Organization, "</Organization>\n", sep="")
cat("      <Pager />\n")
cat("   </User>\n")
cat("   <Location OID=\"LOC.1\" Name=\"", Organization, "\" LocationType=\"Other\">\n", sep="")
cat("      <MetaDataVersionRef StudyOID=\"", StudyOID, "\" MetaDataVersionOID=\"MD.1\" ",
           "EffectiveDate=\"", format(Sys.time(), "%Y-%m-%d"), "\" />\n", sep="") 
cat("   </Location>\r")
cat("</AdminData>\r")
cat("</ODM>\r")

sink()

print("Finished")

}
