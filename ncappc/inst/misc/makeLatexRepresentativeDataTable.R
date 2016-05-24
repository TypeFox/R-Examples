#####################################################################################################################################################################
##                                                                                                                                                                 ##  
## Function for creating latex table over representative portion of data (using threeparttable)                                                                    ##  
## Input: df = "Data frame to display (or NULL if dataDir and dataFile is used)", default=NULL                                                                     ##
##        dataDir   = "The data directory, used to read in data file if df not specified", default=NULL                                                            ##
##        dataFile  = "The data file, used to specify data file if df not specified", default=NULL                                                                 ##
##        exprSubset = "Expression to subset the data set.", default=NULL (no subsetting)                                                                          ##
##        cColumns= "vector of columns to display", default=NULL (all)                                                                                             ##
##        iColWidth = "Maximum number of columns to display on each page", default=15                                                                              ##  
##        iRowNums  = "Maximum number of rows to display", default=20                                                                                              ##  
##        strShortCaption    = "Short caption to use for LOT", default="Representative portion of the data set"                                                    ##
##        strLongCaption   = "Long caption to use for LOT", default="Representative portion of the data set"                                                       ##
##        strFootNote   = "Footnote text to add below the table", default=NULL (no footnote)                                                                       ##
##        strFootNoteOnTables= "If footnotes (strFootNote!=NULL) add them to "all" tables, "first" table, "last" table or "fl" first and last table, default="all" ##
##        cColSeq = "Columns sequence to tabular, default=NULL (automatic).                                                                                        ##
##                   E.g. cColSeq=c(1,9,13,17) will have a table for col 1-9, a table for col 10-13 and a last table for column 13-17                              ##      
##        bShowRowNames = "Show row number/names in table (as leftmost column in each table)", default=TRUE                                                        ##
##        strLabel = "Label of table", default="tab:repData"                                                                                                       ##
##        bReturndf    = "Don't creat any latex output, instead return the subsetted data frame", default=FALSE (produce latex)                                    ##
##                                                                                                                                                                 ##
#####################################################################################################################################################################


makeLatexRepresentativeDataTable <- function(df=NULL,dataDir="",dataFile="",exprSubset=NULL,cColumns=NULL,iColWidth=15,iRowNums=20,
                                             strShortCaption="Representative portion of the data set",
                                             strLongCaption="Representative portion of the data set",
                                             strFootNote=NULL,strFootNoteOnTables="all",cColSeq=NULL,
                                             bShowRowNames=TRUE,strLabel="tab:repData",bReturndf=FALSE){

options(warn = -1)

catt <- function (bPrint,str) if (bPrint) cat(str)  

df <- dataFile

# if (is.null(df)) { #Either read in data set using dataDir and dataFile or using a pre-existing data frame df
#   df <- read.csv(paste(dataDir,dataFile,sep=""))
# }

if (!is.null(df) & dataDir=="" & dataFile=="") {
  print("ERROR: data.frame or dataFile must be specified!!")
}

if (!is.null(exprSubset)) { ## Get a subset of the dataset
  df<-subset(df,eval(exprSubset))
}

#Select some columns only
if (!is.null(cColumns)) {
  df<-df[,cColumns]
}

#Select a few rows
df<-df[1:min(iRowNums,nrow(df)),]


if (is.null(cColSeq)) {

  #Get vector of column sequences
  cColSeq<-seq(0,ncol(df),iColWidth)
  ## If the vector doesn't add up with colwidth
  if (ncol(df)>cColSeq[length(cColSeq)]) cColSeq<-c(cColSeq,ncol(df))
}
else {
  cColSeq[1]<-cColSeq[1]-1
}


for (i in 1:(length(cColSeq)-1)) {
  
  catt(!bReturndf,"\\begin{table}[H]\n")
  catt(!bReturndf,"\\centering\\footnotesize\n")
  catt(!bReturndf,"\\begin{threeparttable}\n")
  
  
  iEnd<-min(ncol(df),cColSeq[i+1])
  df_tmp<-as.data.frame(df[,(cColSeq[i]+1):iEnd,drop = FALSE])
  strColumn<-paste("columns ",cColSeq[i]+1,"-",iEnd,".",sep="")
  if (iEnd==cColSeq[i]+1) strColumn<-paste("column ",iEnd,".",sep="")
  
  
  if (i==1) {
    catt(!bReturndf,paste("\\caption[",strShortCaption,"]{",strLongCaption,", ",strColumn,"}\n",sep=""))
    catt(!bReturndf,paste("\\label{",strLabel,"}\n",sep=""))
  }
  else {
    catt(!bReturndf,paste("\\captionsetup{list=false} \\caption{(continued), ",strColumn,"}\n",sep=""))
  }
  
  
  if (bReturndf==FALSE) print(xtable(df_tmp),include.rownames=bShowRowNames,rotate.colnames=FALSE, size = "footnotesize",floating=FALSE)
  
  if (!is.null(strFootNote) && (strFootNoteOnTables=="all" || 
                                  ((strFootNoteOnTables=="last" || strFootNoteOnTables=="fl") && i==length(cColSeq)-1) ||
                                  ((strFootNoteOnTables=="first" || strFootNoteOnTables=="fl") && i==1))){ ##Add footnotes
    catt(!bReturndf,"\\begin{tablenotes}\n")
    catt(!bReturndf,strFootNote)
    catt(!bReturndf,"\\end{tablenotes}")
  } 
  catt(!bReturndf,"\\end{threeparttable}\n")
  catt(!bReturndf,"\\end{table}\n")
  
  if (i!=length(cColSeq)-1) #Reset the counter except for the last table
  catt(!bReturndf,"\\addtocounter{table}{-1}")
  
}

if (bReturndf==TRUE) return(df)
  
}