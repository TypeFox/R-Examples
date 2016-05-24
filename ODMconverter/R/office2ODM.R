# Martin Dugas 2013
# license: GPL

office2ODM <- function( officefile="") {

if (officefile == "") officefile <- file.choose()

cat("office2ODM (filetypes: csv, txt, xlsx)\n")
cat(paste("officefile =", officefile), "\n")

if ( grepl("\\.xlsx$", officefile ,ignore.case=T) )
{
   csv <- read.xlsx(officefile, 1, header=F, stringsAsFactors=F, colIndex=1:54)
}
else
{   
   if ( grepl("\\.csv$", officefile ,ignore.case=T) )
   {
      csv <- read.csv2(officefile, header=F, as.is=T, colClasses=rep("character",54) )
   }
   else
   {
      if ( grepl("\\.txt$", officefile ,ignore.case=T) )
      {
         csv <- read.delim2(officefile, header=F, as.is=T)
      }
      else stop("unknown file type")
   }
}

# substitute special characters
f1 <- function (x) gsub("&","&amp;",x)
f2 <- function (x) gsub("'","&apos;",x)
f3 <- function (x) gsub("<","&lt;",x)
f4 <- function (x) gsub(">","&gt;",x)
f5 <- function (x) gsub("\"","&quot;",x)
f6 <- function (x) gsub("\x8e","&#196;",x)
f7 <- function (x) gsub("\x99","&#214;",x)
f8 <- function (x) gsub("\x9a","&#220;",x)
f9 <- function (x) gsub("\x84","&#228;",x)
f10 <- function (x) gsub("\x94","&#246;",x)
f11 <- function (x) gsub("\x81","&#252;",x)
f12 <- function (x) gsub("\xe1","&#223;",x)
f13 <- function (x) gsub("\xE4","ae",x)
f14 <- function (x) gsub("\xC2\xB5","&#181;",x)   # micro
f15 <- function (x) gsub("\xFC","ue",x)
csv <- apply(csv,2, f1)
csv <- apply(csv,2, f2)
csv <- apply(csv,2, f3)
csv <- apply(csv,2, f4)
csv <- apply(csv,2, f5)
csv <- apply(csv,2, f6)
csv <- apply(csv,2, f7)
csv <- apply(csv,2, f8)
csv <- apply(csv,2, f9)
csv <- apply(csv,2, f10)
csv <- apply(csv,2, f11)
csv <- apply(csv,2, f12)
csv <- apply(csv,2, f13)
csv <- apply(csv,2, f14)
csv <- apply(csv,2, f15)
csv <- as.data.frame(csv, stringsAsFactors=F)
csv[is.na(csv)] <- ""

StudyOID <- csv[1,2]
Sponsor <- csv[2,2]
Condition <- csv[3,2]
StudyName <- csv[4,2]
StudyDescription <- csv[5,2]
Form <- csv[6,2]
FirstName <- csv[7,2]
LastName <- csv[8,2]
Organization <- csv[9,2]
usedCols <- length(csv[11,])
Context <- csv[11,5:usedCols]

filename <- Form
if (StudyOID != "") filename <- paste(filename, StudyOID, sep="_")
filename <- paste(filename, ".xml", sep="")
# remove \/:*?\"<>| from filename
filename <- gsub("\\","",filename, fixed=T)
filename <- gsub("/","",filename, fixed=T)
filename <- gsub(":","",filename, fixed=T)
filename <- gsub(":","",filename, fixed=T)
filename <- gsub("*","",filename, fixed=T)
filename <- gsub("?","",filename, fixed=T)
filename <- gsub("\"","",filename, fixed=T)
filename <- gsub("<","",filename, fixed=T)
filename <- gsub(">","",filename, fixed=T)
filename <- gsub("|","",filename, fixed=T)
print(paste("Generating", filename))
sink(filename)

cat("<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\r", 
    "<ODM xmlns=\"http://www.cdisc.org/ns/odm/v1.3\"\r", 
    "Description=\"", Form, " ", StudyOID, " ", Condition, "\"\r",
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

cat("<FormDef OID=\"F.1\" Name=\"", Form, "\" Repeating=\"No\">\r",
sep="")
# count itemgroups
n <- dim(csv)[1] - 11   # skip header lines
itemgroup <- ""
aliastags <- ""
IG <- 0
CL <- 0
CLxml <- ""
IGstart <- 1
for (i in 1:n) {
   line <- csv[i+11,]
   line[is.na(line)] <- ""
   line[1] <- tolower(line[1])
   if (line[1] == "itemgroup")
   {   IG <- IG + 1 }
   else if (IG == 0 && nchar(line[1]) >0) { IGstart <- 0 }
}

for (i in IGstart:IG) {
   cat("<ItemGroupRef ItemGroupOID=\"IG.", i, "\" Mandatory=\"Yes\" OrderNumber=\"",i,"\" />\r",
       sep="")
}
cat("</FormDef>\r")

# first items without itemgroup (undefined item group)
if (IGstart == 0) 
{
      cat("<ItemGroupDef OID=\"IG.0\" Name=\"Undefined item group\" Repeating=\"No\">\n",
           "   <Description>\n      <TranslatedText xml:lang=\"en\">Undefined item group</TranslatedText>\n   </Description>\n",
          sep="")
}

IG <- 0
itemno <- 0
for (i in 1:n) {
   line <- csv[i+11,]
   line[is.na(line)] <- ""
   line[1] <- tolower(line[1]) 
   nextline <- csv[i+12,]
   nextline[is.na(nextline)] <- ""
   if (line[1] == "itemgroup") 
   {
      IG <- IG + 1
      if (IG > IGstart) cat(paste(aliastags, "</ItemGroupDef>\n", sep=""))
      aliastags <- ""
      itemno <- 0

      name <- as.character(line[2])
      en <- as.character(line[3])
      if (en == "") en <- name
      de <- as.character(line[4])
      if (de == "") de <- en   
      translated <- paste("      <TranslatedText xml:lang=\"en\">", en, "</TranslatedText>\n",
                          "      <TranslatedText xml:lang=\"de\">", de, "</TranslatedText>\n", sep="")
      cat("<ItemGroupDef OID=\"IG.", IG, "\" Name=\"", name, "\" Repeating=\"No\">\n",
          "   <Description>\n", translated,  "   </Description>\n",
          sep="")
      # Alias for itemgroups 
      Name <- line[5:usedCols]
      for (k in 1:length(Context)) {
         if (Name[k] == "") next
         aliastags <- paste(aliastags, 
                            "   <Alias Context=\"", as.character(Context[k]), 
                            "\" Name=\"", as.character(Name[k]), "\" />\n", sep="" )
      }
      # postcoordination for itemgroups
      while(tolower(nextline[1]) == "postcoordination" ) {
         Name <- nextline[5:usedCols]
         for (k in 1:length(Context)) {
            if (Name[k] == "") next
            aliastags <- paste(aliastags, "   <Alias Context=\"", as.character(Context[k]), 
                "\" Name=\"", as.character(Name[k]), "\" />\n", sep="")
         }
         i <- i + 1
         line <- csv[i+11,]
         line[is.na(line)] <- ""
         line[1] <- tolower(line[1]) 
         nextline <- csv[i+12,]
         nextline[is.na(nextline)] <- ""
      }
   }
   else if (line[1] == "string"  || line[1] == "integer" || line[1] == "float" ||
            line[1] == "boolean" || line[1] == "time"    || line[1] == "date" ) {
      itemno <- itemno + 1
      cat("   <ItemRef ItemOID=\"I.", 1000 * IG + itemno, 
           "\" Mandatory=\"Yes\" OrderNumber=\"", itemno, "\" />\n", sep="")
   }
}
cat(paste(aliastags, "</ItemGroupDef>\n", sep=""))

# ItemDefs
language <- csv[11,4]
if (language == "") { language <- "de" }
IG <- 0
itemno <- 0
for (i in 1:n) 
{
   line <- csv[i+11,] 
   if (line[1] == "itemgroup") 
   {
      IG <- IG + 1
      itemno <- 0
   }
   else if (line[1] == "string"  || line[1] == "integer" || line[1] == "float" ||
            line[1] == "boolean" || line[1] == "time"    || line[1] == "date" ) 
   {
      item <- as.character(line[2])
      itemno <- itemno + 1
      itemlong  <- as.character(line[3])
      itemlong2 <- as.character(line[4])
      if (itemlong  == "") { itemlong   <- item }
      if (itemlong2 == "") { itemlong2 <- itemlong }
      datatype <- as.character(line[1])
      Name <- line[5:usedCols]
      nextline <- csv[i+12,]
      nextline[is.na(nextline)] <- ""
      # Itemdef with OID, Name and DataType
      cat("<ItemDef OID=\"I.", 1000 * IG + itemno, "\" Name=\"", item, 
          "\" DataType=\"", datatype, "\">\r",
          "   <Question>\r",
          "      <TranslatedText xml:lang=\"en\">", itemlong, "</TranslatedText>\r",
          "      <TranslatedText xml:lang=\"", language,"\">", itemlong2, "</TranslatedText>\r",  
          "   </Question>\r", sep="")
      # Alias for items
      currentalias <- ""
      for (k in 1:length(Context)) {
         if (Name[k] == "") next
         currentalias <- paste(currentalias, "   <Alias Context=\"", as.character(Context[k]), 
             "\" Name=\"", as.character(Name[k]), "\" />\n", sep="")
      }
      # postcoordination for items
      while(tolower(nextline[1]) == "postcoordination" ) {
         Name <- nextline[5:usedCols]
         for (k in 1:length(Context)) {
            if (Name[k] == "") next
            currentalias <- paste(currentalias, "   <Alias Context=\"", as.character(Context[k]), 
                "\" Name=\"", as.character(Name[k]), "\" />\n", sep="")
         }
         i <- i + 1
         nextline <- csv[i+12,]
         nextline[is.na(nextline)] <- ""
      }
      # Special CodeList
      if ( tolower(nextline[1]) == "codelist" ) {
         cat("   <CodeListRef CodeListOID=\"", as.character(nextline[2]), "\"/>\n", sep="")
      }
      # Codelistitems
      if ( tolower(nextline[1]) == "codelistitem" ) {
         CL <- CL + 1
         cat("   <CodeListRef CodeListOID=\"CL.", CL, "\"/>\r", sep="")
         CLxml <- paste(CLxml, 
                  "<CodeList OID=\"CL.", CL, "\" Name=\"", item, "\" DataType=\"", datatype, "\">\r",
                  sep="")
         while ( tolower(nextline[1]) == "codelistitem" ) {
            name2 <- as.character(nextline[2])
            en2 <- as.character(nextline[3])
            if (en2 == "") en2 <- name2
            de2 <- as.character(nextline[4])
            if (de2 == "") de2 <- en2   
            CLxml <- paste(CLxml,
                     "   <CodeListItem CodedValue=\"", nextline[2], "\">\n",
                     "      <Decode>\n",
                     "         <TranslatedText xml:lang=\"en\">", en2, "</TranslatedText>\n",
                     "         <TranslatedText xml:lang=\"de\">", de2, "</TranslatedText>\n",
                     "      </Decode>\n", sep="")
            # Alias for codelistitems
            Name <- nextline[5:usedCols]
            for (k in 1:length(Context)) {
               if (Name[k] == "") next
               CLxml <- paste(CLxml, "      <Alias Context=\"", as.character(Context[k]), 
                        "\" Name=\"", as.character(Name[k]), "\" />\n", sep="")
            }
            i <- i + 1
            nextline <- csv[i+12,]
            nextline[is.na(nextline)] <- ""
            # postcoordination for codelistitems
            while(tolower(nextline[1]) == "postcoordination" ) {
               Name <- nextline[5:usedCols]
               for (k in 1:length(Context)) {
                  if (Name[k] == "") next
                  CLxml <- paste(CLxml, "   <Alias Context=\"", as.character(Context[k]), 
                  "\" Name=\"", as.character(Name[k]), "\" />\n", sep="")
               }
               i <- i + 1
               nextline <- csv[i+12,]
               nextline[is.na(nextline)] <- ""
            }
            CLxml <- paste(CLxml, "   </CodeListItem>\n", sep="")
         }
         CLxml <- paste(CLxml, "</CodeList>\n", sep="")
      }
      # output Alias after CodeListRef
      cat(paste(currentalias, "</ItemDef>\n", sep=""))
   }
}

# output Codelists
cat(CLxml)
cat("</MetaDataVersion>\r")
cat("</Study>\r\r")

cat("<AdminData>\r")
cat("   <User OID=\"USR.1\">\r")
cat("      <LoginName>", LastName, "</LoginName>\r", sep="") 
cat("      <DisplayName />\r") 
cat("      <FullName>", FirstName, " ", LastName, "</FullName>\r", sep="") 
cat("      <FirstName>", FirstName, "</FirstName>\r", sep="")
cat("      <LastName>", LastName, "</LastName>\r", sep="") 
cat("      <Organization>", Organization, "</Organization>\r", sep="")
cat("      <Pager />\r")
cat("   </User>\r")
cat("   <Location OID=\"LOC.1\" Name=\"", Organization, "\" LocationType=\"Other\">\r", sep="")
cat("      <MetaDataVersionRef StudyOID=\"", StudyOID, "\" MetaDataVersionOID=\"MD.1\" ",
           "EffectiveDate=\"", format(Sys.time(), "%Y-%m-%d"), "\" />\r", sep="") 
cat("   </Location>\r")
cat("</AdminData>\r")
cat("</ODM>\r")

sink()

cat("Finished\n")

}
