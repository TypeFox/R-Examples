# compareODM
# (C) Martin Dugas 2013
# free for academic use
# input: 2 ODM files
# output: list of identical, matching, similar and differing items

compareODM <- function( ODMfile1="", ODMfile2="") {
   library(XML)
   verbose <- F
   if (ODMfile1 == "") ODMfile1 <- file.choose()
   currentdir <- getwd()
   setwd(dirname(ODMfile1))
   if (ODMfile2 == "") ODMfile2 <- file.choose()

   # parse ODMfile1
   ODM = xmlRoot(xmlTreeParse(ODMfile1, useInternalNodes=TRUE))
   MD <- ODM[["Study"]][["MetaDataVersion"]]
   IDefnodes <- MD[names(xmlChildren(MD)) == "ItemDef"]
   CLnodes <- MD[names(xmlChildren(MD)) == "CodeList"]
   N1 <- vector()
   D1 <- vector()
   U1 <- vector()
   C1 <- vector()
   OID1 <- vector()
   CLOID1 <- vector()

   for (m in 1: length(IDefnodes) ) 
   {
      Inode <- IDefnodes[[m]]
      Name <- xmlAttrs(Inode)["Name"]
      OID <- xmlAttrs(Inode)["OID"]
      DataType <- xmlAttrs(Inode)["DataType"]
      UMLS <- ""
      Alias <- Inode[names(xmlChildren(Inode)) == "Alias"]
      if (length(Alias) > 0) 
	  {
         for (k in 1: length(Alias) ) 
		 {
            Context <- xmlAttrs(Alias[[k]])["Context"]
            if (grepl("^UMLS", Context)) {
               UMLS <- paste(UMLS, xmlAttrs(Alias[[k]])["Name"], " ", sep="")
            }
         }
         UMLS <- sub(" $","",UMLS)
         UMLS <- paste(sort(unlist(strsplit(UMLS, " "))), collapse= " ")
      }
      CodeListRef <- Inode[["CodeListRef"]]
      CodeList <- ""
      found <- F
      if (is.null(CodeListRef)) {   CodeListOID <- ""   }
      else 
	  {
         CodeListOID <- xmlAttrs(CodeListRef)["CodeListOID"]
	     for (k in 1: length(CLnodes) ) 
		 {
            Cnode <- CLnodes[[k]]
            if (CodeListOID == xmlAttrs(Cnode)["OID"]) 
			{
			   found <- T
               CLI <- xmlChildren(Cnode)
               for (l in 1: length(CLI) ) 
			   {
                  CodedValue <- xmlAttrs(CLI[[l]])["CodedValue"]
                  label <- xmlValue(CLI[[l]][["Decode"]][["TranslatedText"]])
                  # determine UMLS codes of codelistitems
                  UMLSCLI <- ""
                  Alias <- CLI[[l]][names(xmlChildren(CLI[[l]])) == "Alias"]
                  if (length(Alias) > 0) {
                     for (x in 1: length(Alias) ) {
                        Context <- xmlAttrs(Alias[[x]])["Context"]
                        if (grepl("^UMLS", Context)) {
                           UMLSCLI <- paste(UMLSCLI, xmlAttrs(Alias[[x]])["Name"], " ", sep="")
                        }
                     }
                  }
                  # sort UMLS codes in codelistitem
                  UMLSCLI <- sub(" $","",UMLSCLI)
                  UMLSCLI <- paste(sort(unlist(strsplit(UMLSCLI, " "))), collapse= " ")
                  CLItext <- paste(CodedValue,"=", label, "(", UMLSCLI, ")", sep="")
                  if (CodeList == "") { CodeList <- CLItext } 
                  else { CodeList <- paste(CodeList, ",", CLItext, sep="") }
               }
            }
         }         
      }
      N1[m] <- Name
      OID1 [m] <- OID
      CLOID1[m] <- CodeListOID
      D1[m] <- DataType
      U1[m] <- UMLS
      C1[m] <- CodeList
	  if (!found) C1[m] <- paste(CodeListOID, C1[m], sep="")
      if (verbose) { cat(Name, DataType, UMLS, CodeList, "\n", sep="|" ) }
   }

   # parse ODMfile2
   ODM = xmlRoot(xmlTreeParse(ODMfile2, useInternalNodes=TRUE))
   MD <- ODM[["Study"]][["MetaDataVersion"]]
   IDefnodes <- MD[names(xmlChildren(MD)) == "ItemDef"]
   CLnodes <- MD[names(xmlChildren(MD)) == "CodeList"]
   N2 <- vector()
   D2 <- vector()
   U2 <- vector()
   C2 <- vector()
   OID2 <- vector()
   CLOID2 <- vector()

   for (m in 1: length(IDefnodes) ) {
      Inode <- IDefnodes[[m]]
      Name <- xmlAttrs(Inode)["Name"]
      OID <- xmlAttrs(Inode)["OID"]
      DataType <- xmlAttrs(Inode)["DataType"]
      UMLS <- ""
      Alias <- Inode[names(xmlChildren(Inode)) == "Alias"]
      if (length(Alias) > 0) {
         for (k in 1: length(Alias) ) {
            Context <- xmlAttrs(Alias[[k]])["Context"]
            if (grepl("^UMLS", Context)) {
               UMLS <- paste(UMLS, xmlAttrs(Alias[[k]])["Name"], " ", sep="")
            }
         }
         # sort UMLS codes
         UMLS <- sub(" $","",UMLS)
         UMLS <- paste(sort(unlist(strsplit(UMLS, " "))), collapse= " ")
      }
      CodeListRef <- Inode[["CodeListRef"]]
      CodeList <- ""
      found <- F
      if (is.null(CodeListRef)) { CodeListOID <- "" }
      else 
	  {
         CodeListOID <- xmlAttrs(CodeListRef)["CodeListOID"]
	     for (k in 1: length(CLnodes) ) 
		 {
            Cnode <- CLnodes[[k]]
            if (CodeListOID == xmlAttrs(Cnode)["OID"]) 
			{
			   found <- T
               CLI <- xmlChildren(Cnode)
               for (l in 1: length(CLI) ) 
			   {
                  CodedValue <- xmlAttrs(CLI[[l]])["CodedValue"]
                  label <- xmlValue(CLI[[l]][["Decode"]][["TranslatedText"]])
                  # determine UMLS codes of codelistitems
                  UMLSCLI <- ""
                  Alias <- CLI[[l]][names(xmlChildren(CLI[[l]])) == "Alias"]
                  if (length(Alias) > 0) 
				  {
                     for (x in 1: length(Alias) ) 
					 {
                        Context <- xmlAttrs(Alias[[x]])["Context"]
                        if (grepl("^UMLS", Context)) 
						{
                           UMLSCLI <- paste(UMLSCLI, xmlAttrs(Alias[[x]])["Name"], " ", sep="")
                        }
                     }
                  }
                  # sort UMLS codes in codelistitem
                  UMLSCLI <- sub(" $","",UMLSCLI)
                  UMLSCLI <- paste(sort(unlist(strsplit(UMLSCLI, " "))), collapse= " ")
                  CLItext <- paste(CodedValue,"=", label, "(", UMLSCLI, ")", sep="")
                  if (CodeList == "") { CodeList <- CLItext } 
                  else { CodeList <- paste(CodeList, ",", CLItext, sep="") }
               }
            }
         }         
      }
      N2[m] <- Name
      OID2[m] <- OID
      CLOID2[m] <- CodeListOID
      D2[m] <- DataType
      U2[m] <- UMLS
      C2[m] <- CodeList
	  if (!found) C2[m] <- paste(CodeListOID, C2[m], sep="")
      if (verbose) { cat(Name, DataType, UMLS, CodeList, "\n", sep="|" ) }
   }


   # compare ODMfile1 and ODMfile2
   cat("ODMfile1=", basename(ODMfile1), 
       "  items: ", length(N1), 
       ", with UMLS code: ", sum(U1 != ""), "\n", sep="")
   cat("ODMfile2=", basename(ODMfile2), 
       "  items: ", length(N2), 
       ", with UMLS code: ", sum(U2 != ""), "\n", sep="")

   identical <- vector()
   similar <- vector()
   matching <- vector()
   status <- rep("d", length(N1))  # default: differing item

   for (i in 1: length(N1) ) 
   {
      if (U1[i] == "") status[i] <- ""   # no UMLS code 

      for (k in 1: length(N2) ) 
	  {
         if (U1[i] == U2[k] && U1[i] != "") 
		 {
            # same concept
            if (D1[i] == D2[k] && C1[i] == C2[k]) 
			{
               # same DataType and same CodeList
               if ( tolower(N1[i]) == tolower(N2[k]) ) 
			   {
                  # same Name
                  identical[i] <- k
                  status[i] <- "i"
               }
               else 
			   {
                  if (U1[i] != "") 
				  {  
                     matching[i] <- k
                     status[i] <- "m"
                  }
               }
            }
            else 
			{ 
               if (U1[i] != "") 
			   {  
                  similar[i] <- k
                  status[i] <- "s"  
               } 
            }
         }
      }
   }

   cat("\nIdentical items: ", sum(!is.na(identical)), "\n", sep="")
   for (i in 1:length(identical) ) {
      if (length(identical) == 0) next
      if (!is.na(identical[i]) ) {
         cat(OID1[i],N1[i],D1[i],U1[i],C1[i],"<->",OID2[identical[i]], "\n", sep="|")
      }
   }

   cat("\nMatching items: ", sum(!is.na(matching)), "\n", sep="")
   for (i in 1:length(matching) ) {
      if (length(matching) == 0) next
      if (!is.na(matching[i]) ) {
         cat(OID1[i],N1[i],D1[i],U1[i],C1[i],"<->",OID2[matching[i]], N2[matching[i]], "\n", sep="|")
      }
   }

   cat("\nSimilar items:", sum(!is.na(similar)), "\n", sep="")
   for (i in 1:length(similar) ) {
      if (length(similar) == 0) next
      if (!is.na(similar[i]) ) {
         cat(U1[i],"\n",sep="")
         cat(OID1[i],N1[i],D1[i],C1[i],"\n", sep="|")
         cat(OID2[similar[i]], N2[similar[i]], D2[similar[i]],C2[similar[i]], "\n", sep="|")
      }
   }

   cat("\nDiffering items:", sum(status == "d"), "\n", sep="")
   for (i in 1:length(status) ) {
      if (status[i] == "d" ) {
         cat(OID1[i],N1[i],D1[i],U1[i],C1[i],"\n", sep="|")
      }
   }

   cat("\nItems without UMLS codes:", sum(status == ""), "\n", sep="")
   for (i in 1:length(status) ) {
      if (status[i] == "" ) {
         cat(OID1[i],N1[i],D1[i],U1[i],C1[i],"\n", sep="|")
      }
   }

   cat("\nFinished\n")
}
