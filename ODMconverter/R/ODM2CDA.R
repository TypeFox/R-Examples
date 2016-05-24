ODM2CDA <- function( ODMfile="", withUMLS=T ) {

if (ODMfile == "") ODMfile <- file.choose()
currentdir <- getwd()
setwd(dirname(ODMfile))

ODM = xmlRoot(xmlTreeParse(ODMfile))
ODMDescription <- xmlAttrs(ODM)["Description"]
CDAfile <- paste(ODMfile,".CDA.xml",sep="")
# cat(paste("CDAfile=",CDAfile,"\r", sep=""))
sink(CDAfile)

StudyName <- xmlValue(ODM[["Study"]][["GlobalVariables"]][["StudyName"]])
StudyDescription <- xmlValue(ODM[["Study"]][["GlobalVariables"]][["StudyDescription"]])
StudyDescription <- gsub("\n"," ",StudyDescription)
StudyDescription <- gsub("\t"," ",StudyDescription)
StudyDescription <- gsub(";",",",StudyDescription)
FirstName <- xmlValue(ODM[["AdminData"]][["User"]][["FirstName"]])
LastName <- xmlValue(ODM[["AdminData"]][["User"]][["LastName"]])
Organization <- xmlValue(ODM[["AdminData"]][["User"]][["Organization"]])

#
# CDA Header
#
cat("<?xml version=\"1.0\"?>
<?xml-stylesheet type=\"text/xsl\" href=\"vhitg-cda-v3.xsl\"?>
<ClinicalDocument xmlns=\"urn:hl7-org:v3\" xmlns:voc=\"urn:hl7-org:v3/voc\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"urn:hl7-org:v3 CDA.xsd\">
   <!-- 
********************************************************
  CDA Header
********************************************************
-->
   <typeId root=\"2.16.840.1.113883.1.3\" extension=\"POCD_HD000040\"/>
   <templateId root=\"2.16.840.1.113883.3.27.1776\"/>
   <id extension=\"c266\" root=\"2.16.840.1.113883.19\"/>
   <code code=\"11488-4\" codeSystem=\"2.16.840.1.113883.6.1\" displayName=\"Consultation note\"/>\r")

cat(paste("   <title>", StudyName, "</title>\r", sep=""))
cat(paste("   <effectiveTime value=\"", format(Sys.time(), "%Y%m%d"), "\"/>\r", sep=""))

cat("   <confidentialityCode code=\"N\" codeSystem=\"2.16.840.1.113883.5.25\"/>
   <setId extension=\"BB35\" root=\"2.16.840.1.113883.19\"/>
   <versionNumber value=\"2\"/>
   <recordTarget>
      <patientRole>
         <id extension=\"12345\" root=\"2.16.840.1.113883.19\"/>
            <patient>
               <name><given>John</given><family>Test</family><suffix>the 2nd</suffix></name>
               <administrativeGenderCode code=\"M\" codeSystem=\"2.16.840.1.113883.5.1\"/>
               <birthTime value=\"19801027\"/>
			</patient>
			<providerOrganization><id extension=\"M345\" root=\"2.16.840.1.113883.19\"/></providerOrganization>
      </patientRole>
   </recordTarget>
   <author>\r")
   
cat(paste("      <time value=\"", format(Sys.time(), "%Y%m%d"), "\"/>\r", sep=""))

cat(paste("      <assignedAuthor>
         <id extension=\"KP00017\" root=\"2.16.840.1.113883.19\"/>
         <assignedPerson>
            <name><given>",FirstName,"</given><family>",LastName,"</family><suffix>MD</suffix></name>
         </assignedPerson>
         <representedOrganization><id extension=\"M345\" root=\"2.16.840.1.113883.19\"/></representedOrganization>
      </assignedAuthor>
   </author>
   <custodian>
      <assignedCustodian>
         <representedCustodianOrganization>
            <id extension=\"M345\" root=\"2.16.840.1.113883.19\"/>
            <name>", Organization ,"</name>
         </representedCustodianOrganization>
      </assignedCustodian>
   </custodian>
   <legalAuthenticator>\r",sep=""))
   
cat(paste("      <time value=\"", format(Sys.time(), "%Y%m%d"), "\"/>\r", sep=""))
   
cat(paste("      <signatureCode code=\"S\"/>
      <assignedEntity>
         <id extension=\"KP00017\" root=\"2.16.840.1.113883.19\"/>
         <assignedPerson><name><given>",FirstName,"</given><family>",LastName,"</family><suffix>MD</suffix></name></assignedPerson>
         <representedOrganization>
            <id extension=\"M345\" root=\"2.16.840.1.113883.19\"/>
         </representedOrganization>
      </assignedEntity>
   </legalAuthenticator>
   <relatedDocument typeCode=\"RPLC\">
      <parentDocument>
         <id extension=\"a123\" root=\"2.16.840.1.113883.19\"/>
         <setId extension=\"BB35\" root=\"2.16.840.1.113883.19\"/>
         <versionNumber value=\"1\"/>
      </parentDocument>
   </relatedDocument>
   <componentOf>
      <encompassingEncounter>
         <id extension=\"KPENC1332\" root=\"2.16.840.1.113883.19\"/>\r",sep=""))
		 
cat(paste("         <effectiveTime value=\"", format(Sys.time(), "%Y%m%d"), "\"/>\r", sep=""))
		 
cat("         <encounterParticipant typeCode=\"CON\">\r")

cat(paste("         <time value=\"", format(Sys.time(), "%Y%m%d"), "\"/>\r", sep=""))

cat(paste("            <assignedEntity>
               <id extension=\"KP00017\" root=\"2.16.840.1.113883.19\"/>
               <assignedPerson><name><given>",FirstName,"</given><family>",LastName,"</family><suffix>MD</suffix></name></assignedPerson>
               <representedOrganization><id extension=\"M345\" root=\"2.16.840.1.113883.19\"/></representedOrganization>
            </assignedEntity>
         </encounterParticipant>
         <location>				
            <healthCareFacility classCode=\"DSDLOC\">
               <code code=\"GIM\" codeSystem=\"2.16.840.1.113883.5.10588\" displayName=\"General internal medicine clinic\"/>
            </healthCareFacility>
         </location>
      </encompassingEncounter>
   </componentOf>\r",sep=""))

#
# CDA Body
#
cat("   <!-- 
********************************************************
  CDA Body
********************************************************
-->
   <component>
      <structuredBody>")

cat("         <!-- 
********************************************************
  Assessment section
********************************************************
-->
         <component>
            <section>
               <code code=\"11496-7\" codeSystem=\"2.16.840.1.113883.6.1\" codeSystemName=\"LOINC\"/>
               <title>Assessment</title>\r")

if (!withUMLS) cat("               <text>\r                  <list>\r")

# output items from ODM form
MD <- ODM[["Study"]][["MetaDataVersion"]]
IGnodes <- MD[names(xmlChildren(MD)) == "ItemGroupDef"]
IDefnodes <- MD[names(xmlChildren(MD)) == "ItemDef"]

newIG <- F
for (i in 1: length(IGnodes) ) {
   itemgroup <- xmlAttrs(IGnodes[[i]])["Name"]
   newIG <- T
   ItemRefNodes <- IGnodes[[i]][names(xmlChildren(IGnodes[[i]]))=="ItemRef"]
   if (length(ItemRefNodes) >0 ) for (k in 1: length(ItemRefNodes)) {
      ItemOID <- xmlAttrs(ItemRefNodes[[k]])["ItemOID"]
      for (m in 1: length(IDefnodes) ) {
         Inode <- IDefnodes[[m]]
         if (ItemOID == xmlAttrs(Inode)["OID"]) {
            itemtext <- xmlAttrs(Inode)["Name"]
            itemtextlong <- itemtext
            tmp <- xmlValue(Inode[["Question"]][["TranslatedText"]])
            if (!is.na(tmp)) { itemtextlong <- tmp }
            DataType <- xmlAttrs(Inode)["DataType"]
            itemgroup <- gsub(";",",",itemgroup)
            if (newIG == T) { newIG <- F }
            itemtext <- gsub(";",",",itemtext)
            itemtextlong <- gsub(";",",",itemtextlong)
            itemtextlong <- gsub("<","&lt;",itemtextlong)
            itemtextlong <- gsub(">","&gt;",itemtextlong)
            if (withUMLS) 
            {
               # with UMLS codes
               cat("               <entry>\r                  <act classCode=\"ACT\" moodCode=\"DEF\">\r")
               UMLScode <- ""
               anodes <- Inode[names(xmlChildren(Inode)) == "Alias"]
               if (length(anodes) > 0 ) for (ii in 1: length(anodes) )
               {
                  anode <- anodes[[ii]]
                  context <- xmlAttrs(anode)["Context"]
                  name <- xmlAttrs(anode)["Name"]
                  if (tolower(substr(context,1,4)) == "umls") {
                     if (UMLScode =="") UMLScode <- name
                     else UMLScode <- paste(UMLScode,name,sep=",")
                  }
               }
               cat(paste("                     <code code=\"", UMLScode, "\" ",
                              "codeSystem=\"2.16.840.1.113883.6.86\" codeSystemName=\"UMLS\" displayName=\"",
                              itemtextlong, "\"/>\r",sep=""))
               cat("                  </act>\r               </entry>\r")
            }
            else {
               # no semantic annotations
               cat(paste("                     <item>", itemtextlong, "</item>\r", sep=""))
            }
          } 
       }	 
   }
}

if (!withUMLS) cat("                  </list>\r               </text>\r")

cat("            </section>
         </component>
      </structuredBody>
   </component>
</ClinicalDocument>\r")


sink()
#cat("Finished\r")
}
