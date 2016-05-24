# Martin Dugas 2013
# license: GPL

ODM2office <- function( ODMfile="", filetype="csv") {

if (ODMfile == "") ODMfile <- file.choose()

ODM = xmlRoot(xmlTreeParse(ODMfile, useInternalNodes=TRUE))
ODMDescription <- xmlAttrs(ODM)["Description"]
StudyOID <- xmlAttrs(ODM[["Study"]])["OID"]
StudyName <- xmlValue(ODM[["Study"]][["GlobalVariables"]][["StudyName"]])
StudyDescription <- xmlValue(ODM[["Study"]][["GlobalVariables"]][["StudyDescription"]])
StudyDescription <- gsub("\n"," ",StudyDescription)
StudyDescription <- gsub("\t"," ",StudyDescription)
StudyDescription <- gsub(";",",",StudyDescription)
Sponsor <- xmlValue(ODM[["Study"]][["GlobalVariables"]][["ProtocolName"]])

MD <- ODM[["Study"]][["MetaDataVersion"]]
Form <- xmlAttrs(MD[["FormDef"]])["Name"]
Condition <- xmlAttrs(MD[["StudyEventDef"]])["Name"]

FirstName <- xmlValue(ODM[["AdminData"]][["User"]][["FirstName"]])
LastName <- xmlValue(ODM[["AdminData"]][["User"]][["LastName"]])
Organization <- xmlValue(ODM[["AdminData"]][["User"]][["Organization"]])

csvfile <- paste(Form, "_", StudyOID,".csv",sep="")
cat("ODMfile=",ODMfile, "\n", sep="")
sink(csvfile)

cat("StudyOID;",StudyOID,"\n", sep="")
cat("Sponsor;",Sponsor,"\n", sep="")
cat("Condition;",Condition,"\n", sep="")
cat("StudyName;",StudyName,"\n", sep="")
cat("StudyDescription;",StudyDescription,"\n", sep="")
cat("Form;",Form,"\n", sep="")
cat("FirstName;",FirstName,"\n",sep="")
cat("LastName;",LastName,"\n",sep="")
cat("Organization;",Organization,"\n;\n",sep="")

IGnodes <- MD[names(xmlChildren(MD)) == "ItemGroupDef"]
IDefnodes <- MD[names(xmlChildren(MD)) == "ItemDef"]
CLnodes <- MD[names(xmlChildren(MD)) == "CodeList"]
Contexts <- c("Type","Name","en","de")
output <- ""

newIG <- F
for (i in 1: length(IGnodes) )
{
   itemgroup <- xmlAttrs(IGnodes[[i]])["Name"]
   itemgroup <- gsub(";",",",itemgroup)
   itemgrouplong <- ""
   itemgroup_de <- ""
   tmp <- xmlValue(IGnodes[[i]][["Description"]][["TranslatedText"]])
   if (!is.na(tmp))
   {
      itemgrouplong <- tmp
      if ( length(xmlChildren(IGnodes[[i]][["Description"]]))> 1 ) itemgroup_de <- xmlValue(IGnodes[[i]][["Description"]][2][["TranslatedText"]])
   }
   output_itemgroup <- c("itemgroup", itemgroup, itemgrouplong, itemgroup_de)

   # alias itemgroup
   anodes <- IGnodes[[i]][names(xmlChildren(IGnodes[[i]])) == "Alias"]
   if (length(anodes) > 0 ) for (ii in 1: length(anodes) )
   {
      anode <- anodes[[ii]]
      context <- xmlAttrs(anode)["Context"]
      name <- xmlAttrs(anode)["Name"]
      if (sum(Contexts == context) == 0) Contexts <- c(Contexts, context) 
      if ( is.na(output_itemgroup[Contexts == context]) ) { output_itemgroup[Contexts == context] <- name }
      else { output_itemgroup[Contexts == context] <- paste(output_itemgroup[Contexts == context], name, sep="   ") }
   } 
   output_itemgroup[is.na(output_itemgroup)] <- ""
   output_itemgroup <- paste(output_itemgroup, collapse=";", sep="")  
   
   newIG <- T
   ItemRefNodes <- IGnodes[[i]][names(xmlChildren(IGnodes[[i]]))=="ItemRef"]
   ### avoid error with empty itemgroups
   if (length(ItemRefNodes) > 0) for (k in 1: length(ItemRefNodes))
   {
      ItemOID <- xmlAttrs(ItemRefNodes[[k]])["ItemOID"]
      for (m in 1: length(IDefnodes) )
      {
         Inode <- IDefnodes[[m]]
         if (ItemOID == xmlAttrs(Inode)["OID"])
         {
            itemtext <- xmlAttrs(Inode)["Name"]
            itemtextlong <- ""
            itemtext_de <- ""
            tmp <- xmlValue(Inode[["Question"]][["TranslatedText"]])
            if (!is.na(tmp)) {
               itemtextlong <- tmp
               if ( length(xmlChildren(Inode[["Question"]]))> 1 ) {
                  itemtext_de <- xmlValue(Inode[["Question"]][2][["TranslatedText"]])
               }
            }
            DataType <- xmlAttrs(Inode)["DataType"]
            if (newIG == T) { newIG <- F; output <- paste(output, output_itemgroup, "\n", sep="") }
            itemtext <- gsub(";",",",itemtext)
            itemtextlong <- gsub(";",",",itemtextlong)
            itemtext_de <- gsub(";",",",itemtext_de)
            tmp <- c(DataType, itemtext, itemtextlong, itemtext_de)
            # alias items
            anodes <- Inode[names(xmlChildren(Inode)) == "Alias"]
            if (length(anodes) > 0 ) for (ii in 1: length(anodes) )
            {
               anode <- anodes[[ii]]
               context <- xmlAttrs(anode)["Context"]
               name <- xmlAttrs(anode)["Name"]
               if (sum(Contexts == context) == 0) Contexts <- c(Contexts, context) 
               if ( is.na(tmp[Contexts == context]) ) { tmp[Contexts == context] <- name }
               else
               { 
                  tmp[is.na(tmp)] <- ""
                  tmp2 <- paste(tmp, collapse=";", sep="") 
                  output <- paste(output, tmp2,"\n", sep="")
                  tmp <- c("postcoordination","","","")
                  tmp[Contexts == context] <- name 
               }
            }
            tmp[is.na(tmp)] <- ""
            tmp2 <- paste(tmp, collapse=";", sep="") 
            output <- paste(output, tmp2,"\n", sep="")

            # codelistitems
            if (!is.null(Inode[["CodeListRef"]])) {
               CLOID <- xmlAttrs(Inode[["CodeListRef"]])["CodeListOID"]
               found <- F
               for (ii in 1: length(CLnodes) )
               {
                  CLnode <- CLnodes[[ii]]
                  if (CLOID == xmlAttrs(CLnode)["OID"])
                  {
                     found <- T
                     CLInodes <- CLnode[names(xmlChildren(CLnode)) == "CodeListItem"]
                     for (iii in 1:length(CLInodes) )
                     {
                        CLInode <- CLInodes[[iii]]
                        CodedValue <- xmlAttrs(CLInode)["CodedValue"]
                        text1 <- ""
                        text2 <- ""
                        tmp <- xmlValue(CLInode[["Decode"]][["TranslatedText"]])
                        if (!is.na(tmp))
                        { 
                           text1 <- tmp
                           if ( length(xmlChildren(CLInode[["Decode"]]))> 1 ) 
                           {
                              text2 <- xmlValue(CLInode[["Decode"]][2][["TranslatedText"]])
                           }
                        }
				# Alias for codelistitems
                        tmp <- c("codelistitem", CodedValue, text1,text2)
                        anodes <- CLInode[names(xmlChildren(CLInode)) == "Alias"]
                        if (length(anodes) > 0 ) for (ii in 1: length(anodes) )
                        {
                           anode <- anodes[[ii]]
                           context <- xmlAttrs(anode)["Context"]
                           name <- xmlAttrs(anode)["Name"]
                           if (sum(Contexts == context) == 0) Contexts <- c(Contexts, context) 

                           if ( is.na(tmp[Contexts == context]) ) { tmp[Contexts == context] <- name }
                           else { tmp[Contexts == context] <- paste(tmp[Contexts == context], name, sep="   ") }
                        }
                        tmp[is.na(tmp)] <- ""
                        tmp <- paste(tmp, collapse=";", sep="") 
                        output <- paste(output, tmp, "\n", sep="")
                     }
                  }
               }
               if (!found) output <- paste(output, "codelist;",CLOID,"\n",sep="")
            }
         }		 
      }	 
   }
}

cat(paste(Contexts,collapse=";",sep=""))
cat("\n")
cat(output)
sink()

if (filetype == "csv")   cat("Output:",csvfile,"\n")

if (filetype == "xlsx") 
{
   suppressWarnings( csv <- read.csv2(csvfile, header=F, as.is=T,colClasses=rep("character",54) ) )
   xlsxfile <- paste(Form, "_", StudyOID, ".xlsx",sep="")
   cat("Output:",xlsxfile,"\n")
   write.xlsx(csv, xlsxfile, row.names=F, col.names=F)
   unlink(csvfile)
}

if (filetype == "txt") 
{
   suppressWarnings( csv <- read.csv2(csvfile, header=F, as.is=T,colClasses=rep("character",54) ) )
   txtfile <- paste(Form, "_", StudyOID,".txt",sep="")
   cat("Output:",txtfile,"\n")
   write.table(csv, file=txtfile, sep="\t", quote=F, row.names=F, col.names=F)
   unlink(csvfile)
}


cat("Finished\n")

}
