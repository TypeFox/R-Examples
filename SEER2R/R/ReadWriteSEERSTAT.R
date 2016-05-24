
#' SEER2R is a package used for reading data from and writing data to SEET*STAT data format
#' SEER*STAT plain data are represented in two files: dic and txt files. The dic file contains 
#'    variable information and the txt file contains data in columns.
NA;
#' This package is used for 
#' @rdname SEER2-package
#roxygen();


######################
# private functions
######################
#' accessory functions
.string.trimleadingtrailing <- function(string){
#    input   see the parameter list
#    return  result

  result <-  sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", string, perl=TRUE)

}

.string.cutSubStringFileExtension <- function(FileName){
#    input   see the parameter list
#    return  result
  ExtRegularExpr <- c("[.][^.]*$");
  result <-  sub(ExtRegularExpr, "", FileName, perl=TRUE);
  result;
}

.catStr=function(myStr, ...,sep=""){
  result=paste(
    paste(myStr,collapse="",sep=sep),
    ...,
    sep=sep  
  );
}

.create.emptyDataFrame.string=function(ColNames){

  TableContentData=data.frame();
  NumCols=length(ColNames);
  if(NumCols>0)
  for(ColIndex in 1:NumCols){
    TableContentData=cbind(TableContentData,data.frame(numeric(0)));
  }
  colnames(TableContentData)=ColNames;
  TableContentData;
}
.is.empty=function(x){
  if(class(x) == "data.frame"){
    result=(dim(x)[1]*dim(x)[2] == 0)
  }else{
    result=(length(as.vector(x)) == 0);
  }
  result;
}

.removeACharAll = function(char, string){
  result <-  gsub(char, "\\1", as.character(string), perl=TRUE)
}

.removeQuoteAll = function(string){
  result <-  .removeACharAll("\"",string);
}

.addQuoteIfSpaceAtEnd = function(string){
  string = as.character(string);
  Pattern1 = "^[[:space:]]+";
  Pattern2 = "[[:space:]]+$";
  ReturnedStr = string;
  SpaceAtEnd = (length(grep(Pattern1, string)) == 1) || (length(grep(Pattern2, string)) == 1);
  if(SpaceAtEnd){
    ReturnedStr = .catStr("\"",string,"\"");
  }
  return(ReturnedStr);
}

.is.substring = function(x, string){
  return(length(grep(x, string,ignore.case = TRUE)) == 1);
}

.substituteColNames = function(ByVarNamesInF, ByVarsDataColNamesInF, ColNamesInF){
  NumNamesInF = length(ColNamesInF);
  NumByVarNamesInF = length(ByVarNamesInF);
  ColIndexByVar = numeric(NumByVarNamesInF);
  if(NumNamesInF > 0 && NumByVarNamesInF > 0)
  for(VInd in 1:NumByVarNamesInF){
    for(NameInd in 1:NumNamesInF){
#      if(.is.substring(ByVarNamesInF[VInd], gsub("_"," ", ColNamesInF[NameInd]))){
      if(.is.substring(ByVarNamesInF[VInd], ColNamesInF[NameInd])){
        ColNamesInF[NameInd] = ByVarsDataColNamesInF[VInd];
        ColIndexByVar[VInd] = NameInd;
        break;
      }
    }
  }
  return(list(ColNamesInF = ColNamesInF,ColIndexByVar = ColIndexByVar));
}

.assignColNames = function(data,ColNames = NULL){
  OldColNames = colnames(data);
  if(!is.null(ColNames)){
    NewColNameInfo = .substituteColNames(ColNames, ColNames, OldColNames);
    colnames(data) = NewColNameInfo$ColNamesInF;
    attr(data,"NewColNameInfo") = NewColNameInfo;
  }
  return(data);
}

.getSubDataByVarName = function(data,ColNames = NULL){
  TempData = .assignColNames(data,ColNames = ColNames);
  NewColNameInfo = attr(TempData,"NewColNameInfo");
  GoodColIndics =  NewColNameInfo$ColIndexByVar[NewColNameInfo$ColIndexByVar>0];
  NumByVars = length(GoodColIndics);
  if(NumByVars > 0){
    TempData1 = TempData[,GoodColIndics];
    #attr(TempData1,"DICInfo") = attr(TempData,"DICInfo");
  }else{
    TempData1 = TempData;
  }
  return(TempData1);
}

#############################
# package specific functions
#############################

.replaceLabOrNum = function(VarFormatSecListInF, dataInF, LabReplaced = TRUE){
#LabReplaced = TRUE: means that current values are Labels, and will be replaced with numeric values
#columns of  VarFormatSecListInF[[VarIndex]]: Numeric, Labels
  NumVars = length(VarFormatSecListInF);
  if(NumVars > 0){
    for(VarIndex in 1:NumVars){
      cntSection = VarFormatSecListInF[[VarIndex]];
      if(!.is.empty(cntSection)){
        NumLabels = dim(cntSection)[1];
        for(LabIndex in 1:NumLabels){
          ValueOnCheck = as.character(cntSection[LabIndex,1]);
          ValueAssigned = .removeQuoteAll(cntSection[LabIndex,2]);
          if(LabReplaced){
            ValueOnCheck = .removeQuoteAll(as.character(cntSection[LabIndex,2]));
            ValueAssigned = cntSection[LabIndex,1];
          }
          dataInF[as.character(dataInF[,VarIndex]) == ValueOnCheck,VarIndex] = ValueAssigned; 
        }   
      }
    }
  }
  return(dataInF);
}

#' @include myLib.r
#' @callGraph
.read.SeerStatDIC = function(DICfileName){
#   input: 
#      DICfileName: 
#   returned: (see the end of the function definition)

  DicFileInF = DICfileName;
	lines = readLines(DicFileInF);
	
	#function for reading a section
	findSectionInfo = function(lines,Pattern,splitchar="=",addedToName="",IsPatternStr=FALSE) {
	#if IsPatternStr=FALSE, Pattern is regular expression
    ItemInfo=data.frame("ItemNameInDic"=numeric(0),"ItemValueInDic"=numeric(0));
    SectionName=character();
    if(length(lines)>0)
  	for(LineInd in 1:length(lines)){
  		if(length(grep(Pattern, lines[LineInd], fixed=IsPatternStr)) == 1){
  			ItemInd = 1;
  			ItemLineInd=1;
  			SectionName=.string.trimleadingtrailing(lines[LineInd]);
  			while(lines[LineInd + ItemLineInd] != "" ){
          if(length(grep(splitchar, lines[LineInd + ItemLineInd])) == 1){
    				ItemNames = unlist(strsplit(lines[LineInd + ItemLineInd], split=splitchar));
    				if(length(ItemNames) == 2){
    					ItemInfo[ItemInd,] = c(paste(as.character(ItemNames[1]),addedToName,sep=""),as.character(ItemNames[2]));
    				}else{
    					ItemInfo[ItemInd,1] = paste(as.character(ItemNames[1]),addedToName,sep="");
    				}
     				ItemInd = ItemInd + 1;
    			}
    			ItemLineInd = ItemLineInd + 1;
  			}
  		}
  	}
  	attr(ItemInfo, "SectionNameInDIC") = SectionName;
  	return(ItemInfo);
	}
	
	#read sections:
  #   System 
	SystemInfo = findSectionInfo(lines,Pattern="System[]]",splitchar="=",addedToName="");
  #   Session Options 
	SessionOptionInfo = findSectionInfo(lines,Pattern="Session Options[]]",splitchar="=",addedToName="");
  #   Export Options 
	ExportOptionInfo = findSectionInfo(lines,Pattern="Export Options[]]",splitchar="=",addedToName="");
  #   VarAllInfo 
	VarAllInfo = findSectionInfo(lines,Pattern="Variables[]]",splitchar="=",addedToName="");

  #   Variables of non-survival data, Page Variables of survival data
	VarInfo = findSectionInfo(lines,Pattern="Variables[]]",splitchar="Name=",addedToName="Name");
  VarBaseInfo = findSectionInfo(lines,Pattern="Variables[]]",splitchar="Base=",addedToName="Base");
  VarFormatSecList = list();
  VarDimInfo = findSectionInfo(lines,Pattern="Variables[]]",splitchar="Dimension=",addedToName="Dimension");
  VarDimList = list();

  findVarFormatSecList = function(lines,VarBaseInfoInF,VarInfoInF){
    NumBases = dim(VarBaseInfoInF)[1];
    PageVar="Page type";
    NumVars = dim(VarInfoInF)[1];
    VarFormatSecListInF = list();

    if(NumVars > 0)
    for(VarInd in NumVars:1){
      cntVar = as.character(VarInfoInF[VarInd,"ItemNameInDic"]);
      cntVarItemValueInDic = as.character(VarInfoInF[VarInd,"ItemValueInDic"]);

      cntPattern = cntVarItemValueInDic;
      cntPattern = paste(cntPattern,"]",sep="");
      VarFormatSecListInF[[VarInd]] = findSectionInfo(lines,Pattern=cntPattern,splitchar="=",IsPatternStr=T);
      if(.is.empty(VarFormatSecListInF[[VarInd]]))
        VarFormatSecListInF[[VarInd]] = NULL;
      comments = function(){
      #if the variable is Page type, then handle it
      if(length(grep(PageVar, cntVarItemValueInDic)) == 1){
        cntPattern = cntVarItemValueInDic;
        cntPattern = paste(cntPattern,"]",sep="");
        VarFormatSecListInF[[VarInd]] = findSectionInfo(lines,Pattern=cntPattern,splitchar="=",IsPatternStr=T);
      }
      if(NumBases > 0)
      for(BaseInd in 1:NumBases){
        cntBase = as.character(VarBaseInfoInF[BaseInd,"ItemNameInDic"]);
        cntBase = unlist(strsplit(cntBase,"Base"))[1];
        if(length(grep(cntBase, cntVar)) == 1){
          cntPattern = cntVarItemValueInDic;
          cntPattern = paste(cntPattern,"]",sep="");
          VarFormatSecListInF[[VarInd]] = findSectionInfo(lines,Pattern=cntPattern,splitchar="=",IsPatternStr=T);
        }
      }
      }
    }
    return(VarFormatSecListInF)
  }
  VarFormatSecList = findVarFormatSecList(lines,VarBaseInfo,VarInfo);

	PageVarNameInReadData = gsub("[,:()<>={}!@#$%^&*+-]", "", VarInfo[,"ItemValueInDic"]);
	PageVarNameInReadData = gsub("[[:space:]]+", "_", PageVarNameInReadData);
	Vars = data.frame(
    "VarSymbolInDic" = VarInfo[,"ItemNameInDic"],
    "VarNameInReadData" = PageVarNameInReadData,
    "VarNameInTXT" = VarInfo[,"ItemValueInDic"]
  );

  attr(Vars, "SectionNameInDIC") = attr(VarInfo, "SectionNameInDIC");
  attr(Vars, "VarBaseInfo") = VarBaseInfo;
  attr(Vars, "VarFormatSecList") = VarFormatSecList;
  attr(Vars, "VarDimInfo") = VarDimInfo;
  attr(Vars, "VarDimList") = VarDimList;

	returned=list(
    SystemInfo=SystemInfo,
    SessionOptionInfo = SessionOptionInfo,
    ExportOptionInfo = ExportOptionInfo,
    VarAllInfo = VarAllInfo,
    
    Vars = Vars,
    VarBaseInfo = VarBaseInfo,
    VarFormatSecList = VarFormatSecList,
    VarDimInfo = VarDimInfo,
    VarDimList = VarDimList
  );
}

####
##
.getReadingOptionFromExportOption = function(ExportOptionInfoInF){
#   header = FALSE, sep = "", quote = "\"'",
#   na.strings = "NA"
#   gzfile("file.text.gz")
  rownames(ExportOptionInfoInF)=ExportOptionInfoInF$ItemNameInDic;
  GZipped = FALSE;
  if("GZipped" %in% ExportOptionInfoInF$ItemNameInDic){
    GZipped = ifelse(ExportOptionInfoInF["GZipped",2] == "true", TRUE, FALSE);
  }

  VariableformatNum = FALSE;
  quotedlabels = TRUE;
  if("Variable format" %in% ExportOptionInfoInF$ItemNameInDic){
    VariableformatNum = !(length(grep("label", ExportOptionInfoInF["Variable format",2])) == 1);
    quotedlabels = ifelse(ExportOptionInfoInF["Variable format",2] == "quotedlabels", TRUE, FALSE);
  }

  FileformatIsDOS = TRUE;
  if("File format" %in% ExportOptionInfoInF$ItemNameInDic){
    FileformatIsDOS = ifelse(ExportOptionInfoInF["File format",2] == "DOS/Windows", TRUE, FALSE);
  }

  Fielddelimiter = "\t";
  if("Field delimiter" %in% ExportOptionInfoInF$ItemNameInDic){
    Fielddelimiter = ifelse(ExportOptionInfoInF["Field delimiter",2] == "space", " ", 
      ifelse(ExportOptionInfoInF["Field delimiter",2] == "comma", ",",
        ifelse(ExportOptionInfoInF["Field delimiter",2] == "semi-colon", ";", "\t")        
      ));
  }

  MissingCh = ".";
  if("Missing character" %in% ExportOptionInfoInF$ItemNameInDic){
    MissingCh = ifelse(ExportOptionInfoInF["Missing character",2] == "period", ".", " ");
  }

  FieldsValuewithdelimiterinquotes = FALSE;
  if("Fields with delimiter in quotes" %in% ExportOptionInfoInF$ItemNameInDic){
    FieldsValuewithdelimiterinquotes = ifelse(ExportOptionInfoInF["Fields with delimiter in quotes",2] == "true", TRUE, FALSE);
  }

  VarNameIncluded = TRUE;
  if("Variable names included" %in% ExportOptionInfoInF$ItemNameInDic){
    VarNameIncluded = ifelse(ExportOptionInfoInF["Variable names included",2] == "true", TRUE, FALSE);
  }
  
  return(list(
    GZipped = GZipped,
    VariableformatNum = VariableformatNum,
    quotedlabels = quotedlabels,
    FileformatIsDOS = FileformatIsDOS,
    Fielddelimiter = Fielddelimiter,
    MissingCh = MissingCh,
    FieldsValuewithdelimiterinquotes = FieldsValuewithdelimiterinquotes,
    VarNameIncluded = VarNameIncluded
  ));

}



#' @include myLib.r
#' @callGraph
.read.SeerStatTXT = function(DICfileName, TXTfileName, UseVarLabelsInData=FALSE,ReadHeaderOnly=FALSE,...){
#   input: 
#      DICfileName, TXTfileName 
#   returned: (see the end of the function definition)

  InfoFromDIC=.read.SeerStatDIC(DICfileName);
  ReadingOptionInfo = .getReadingOptionFromExportOption(InfoFromDIC$ExportOptionInfo);
  NumVars = length(InfoFromDIC$VarFormatSecList);
  if(NumVars == 0)
    UseVarLabelsInData = TRUE;

  DICInfo = list(
    SystemInfo = InfoFromDIC$SystemInfo,
    SessionOptionInfo = InfoFromDIC$SessionOptionInfo,
    ExportOptionInfo = InfoFromDIC$ExportOptionInfo,
    VarAllInfo = InfoFromDIC$VarAllInfo,
    VarFormatSecList = InfoFromDIC$VarFormatSecList,
    
    UseVarLabelsInData = UseVarLabelsInData
  )
  #DataFile
  if(!ReadHeaderOnly){
    if(ReadingOptionInfo$GZipped){
      if(is.null(TXTfileName)){
        TXTfileName = paste(.string.cutSubStringFileExtension(DICfileName), ".gz",sep="");
      }
      DataFile = gzfile(TXTfileName);  
    }else{
      if(is.null(TXTfileName)){
        TXTfileName = paste(.string.cutSubStringFileExtension(DICfileName), ".txt",sep="");
      }
      DataFile = TXTfileName;  
    }
    
    #cntquote
    cntquote = ifelse(!ReadingOptionInfo$VariableformatNum && ReadingOptionInfo$quotedlabels || ReadingOptionInfo$VarNameIncluded, "\"", "");  
    NumSkippedLines = ifelse(ReadingOptionInfo$VarNameIncluded,1,0);
    TXTData = read.table(DataFile, skip = NumSkippedLines, #header = ReadingOptionInfo$VarNameIncluded,
      sep = ReadingOptionInfo$Fielddelimiter, na.strings = ReadingOptionInfo$MissingCh,
      quote = cntquote, ...   
      );
      
    colnames(TXTData) = InfoFromDIC$Vars$VarNameInReadData;
    if(UseVarLabelsInData && NumVars > 0){
      TXTData = .replaceLabOrNum(InfoFromDIC$VarFormatSecList, TXTData, LabReplaced = FALSE);
    }
  }else{
    return(DICInfo);
  }
  attr(TXTData, "DICInfo") = DICInfo;
  attr(TXTData, "assignColNames") = .assignColNames;
  attr(TXTData, "getSubDataByVarName") = .getSubDataByVarName;
  return(TXTData);
}


.write.SeerStatDIC.section = function(DICSection, filename, append=TRUE){
#Assume that DICSection has very thing required ready
  SectionNameInDIC = attr(DICSection, "SectionNameInDIC");
  outputStr = .catStr(SectionNameInDIC,"\n");
  
  NumRows = dim(DICSection)[1];
  if(NumRows > 0)
  for(RowIndex in 1:NumRows){
    outputStr = .catStr(outputStr,
      DICSection$ItemNameInDic[RowIndex],
      "=",
      DICSection$ItemValueInDic[RowIndex],
      "\n"
    );
  }
  write(outputStr, file=filename, sep="",append=append);
}

.assignItemValueToSecInfoVar = function(cntItemValue,cntItemName, cntSection){
  SectionColNames = c("ItemNameInDic","ItemValueInDic");
  if(cntItemName %in% cntSection$ItemNameInDic){
    cntSection[cntItemName == cntSection$ItemNameInDic,2] = cntItemValue;
  }else{
    cntRow = .create.emptyDataFrame.string(SectionColNames);
    cntRow[1,1] = cntItemName;
    cntRow[1,2] = cntItemValue;

    cntSectionNameInDIC = attr(cntSection, "SectionNameInDIC");
    cntSection = rbind(cntRow, cntSection);
    attr(cntSection, "SectionNameInDIC") = cntSectionNameInDIC; 
  }
  return(cntSection);
}

.write.SeerStatDIC = function(DICInfo, DICfileName){
#   input: 
#      DICfileName, TXTfileName 
#   returned: 
#     a list storing info in DICfileName
  .write.SeerStatDIC.section(DICInfo$SystemInfo, DICfileName, append=FALSE)
  .write.SeerStatDIC.section(DICInfo$SessionOptionInfo, DICfileName, append=TRUE)
  .write.SeerStatDIC.section(DICInfo$ExportOptionInfo, DICfileName, append=TRUE)
  .write.SeerStatDIC.section(DICInfo$VarAllInfo, DICfileName, append=TRUE)
  
  NumMaxFormatSec = length(DICInfo$VarFormatSecList);
  if(NumMaxFormatSec > 0)
  for(VarCaseInd in 1:NumMaxFormatSec){
    if(!.is.empty(DICInfo$VarFormatSecList[[VarCaseInd]]))
      .write.SeerStatDIC.section(DICInfo$VarFormatSecList[[VarCaseInd]], DICfileName, append=TRUE);
  }

  return(DICInfo);
}

.write.SeerStatTXT = function(myData, DICInfo,...){

  ReadingOptionInfo = .getReadingOptionFromExportOption(DICInfo$ExportOptionInfo);

  #DataFile
  TXTfileName = DICInfo$SystemInfo["Output filename" == DICInfo$SystemInfo$ItemNameInDic,2];
  if(ReadingOptionInfo$GZipped){
    DataFile = gzfile(TXTfileName);  
  }else{
    DataFile = TXTfileName;  
  }
  #cntquote
  cntquote = (!ReadingOptionInfo$VariableformatNum && ReadingOptionInfo$quotedlabels);
  outputColNames = colnames(myData);
  if(!cntquote)
    outputColNames = paste("\"",outputColNames,"\"",sep="");
  if(!ReadingOptionInfo$VarNameIncluded){
    outputColNames = FALSE;
  }
  
  write.table(myData, DataFile, col.names = outputColNames,
    sep = ReadingOptionInfo$Fielddelimiter, na = ReadingOptionInfo$MissingCh,
    quote = cntquote,
    row.names = FALSE,...   
    );
}

.prepareDICInfoAndData = function(dataInF,DICfileName,TXTfileName=NULL,UseVarLabelsInTxtFile=TRUE,LabVarsNames=NULL){
  
  DICInfo = attr(dataInF, "DICInfo");
  SectionColNames = c("ItemNameInDic","ItemValueInDic");
  if(is.null(DICInfo)){
    DICInfo = list();
    DICInfo$SystemInfo = .create.emptyDataFrame.string(SectionColNames);
    attr(DICInfo$SystemInfo, "SectionNameInDIC") = "[System]";
    DICInfo$SessionOptionInfo = .create.emptyDataFrame.string(SectionColNames);
    attr(DICInfo$SessionOptionInfo, "SectionNameInDIC") = "[Session Options]";
    DICInfo$ExportOptionInfo = .create.emptyDataFrame.string(SectionColNames);
    attr(DICInfo$ExportOptionInfo, "SectionNameInDIC") = "[Export Options]";
    DICInfo$VarAllInfo = .create.emptyDataFrame.string(SectionColNames);
    attr(DICInfo$VarAllInfo, "SectionNameInDIC") = "[Life Page Variables]";
    checkSurvData = function(ColNames){
      IsSurvDataInF = FALSE;
      NumNames = length(ColNames);
      if(NumNames > 0)
      for(NameIndex in 1:NumNames){
        cntName = ColNames[NameIndex];
        if( (length(grep("Page.*type", cntName)) == 1) ){
          IsSurvDataInF = TRUE;
          break;            
        }
      }
      return(IsSurvDataInF);    
    }
    IsSurvData = checkSurvData(colnames(dataInF));
    if( !IsSurvData )
      attr(DICInfo$VarAllInfo, "SectionNameInDIC") = "[Variables]";
        
    DICInfo$VarFormatSecList = list();
  
    if(is.null(TXTfileName)){
      TXTfileName = paste(.string.cutSubStringFileExtension(DICfileName), ".txt",sep="");
    }
    DICInfo$SystemInfo[1,1] = "Output filename";
    DICInfo$SystemInfo[1,2] = TXTfileName;
    

    #DICInfo$ExportOptionInfo
    RowInd = 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "GZipped";
    DICInfo$ExportOptionInfo[RowInd,2] = "false";
    
    RowInd = RowInd + 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "Variable format";
    DICInfo$ExportOptionInfo[RowInd,2] = "quotedlabels";
    if(!UseVarLabelsInTxtFile)
      DICInfo$ExportOptionInfo[RowInd,2] = "numeric";
    
    RowInd = RowInd + 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "Field delimiter";
    DICInfo$ExportOptionInfo[RowInd,2] = "tab";
    
    RowInd = RowInd + 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "Missing character";
    DICInfo$ExportOptionInfo[RowInd,2] = "period";
    
    RowInd = RowInd + 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "Fields with delimiter in quotes";
    DICInfo$ExportOptionInfo[RowInd,2] = "true";
    
    RowInd = RowInd + 1;
    DICInfo$ExportOptionInfo[RowInd,1] = "Variable names included";
    DICInfo$ExportOptionInfo[RowInd,2] = "false";
    
    #DICInfo$VarAllInfo
    NumVars = dim(dataInF)[2];
    VarNames = colnames(dataInF);
    if(NumVars > 0)
    for(VarIndex in 1:NumVars){
      DICInfo$VarAllInfo[VarIndex,1] = paste("Var",VarIndex,"Name",sep="");
      DICInfo$VarAllInfo[VarIndex,2] = VarNames[VarIndex];
    }
    if(!UseVarLabelsInTxtFile){
      NumVars = dim(dataInF)[2];
      VarNames = colnames(dataInF);
      if(!is.null(LabVarsNames)){
        VarsNamesOnCheck = LabVarsNames;
      }else{
        VarsNamesOnCheck = VarNames;
      }
      NumLabVars = length(VarsNamesOnCheck);
      if(NumVars > 0 && NumLabVars > 0)
      for(VarIndex in 1:NumVars){
        for(LabVIndex in 1:NumLabVars){
          if(as.character(VarsNamesOnCheck[LabVIndex]) == as.character(VarNames[VarIndex]) && 
            (class(dataInF[,VarIndex]) == "character" || class(dataInF[,VarIndex]) == "factor")){
            DistinctValues = dataInF[,VarIndex];
            DistinctValues = levels(factor(DistinctValues));
            NumValues = length(DistinctValues);
            if(NumValues > 0){
              DICInfo$VarFormatSecList[[VarIndex]] = .create.emptyDataFrame.string(SectionColNames);
              for(ValueIndex in 1:NumValues){
                DICInfo$VarFormatSecList[[VarIndex]][ValueIndex,1] = ValueIndex -1;
                DICInfo$VarFormatSecList[[VarIndex]][ValueIndex,2] = .addQuoteIfSpaceAtEnd(DistinctValues[ValueIndex]);
              }
              cntSectionName = paste("[Format=",VarsNamesOnCheck[LabVIndex],"]",sep="");
              attr(DICInfo$VarFormatSecList[[VarIndex]], "SectionNameInDIC") = cntSectionName;
            }
          }
        }
      }
    }   

  }else{
    ReadingOptionInfo = .getReadingOptionFromExportOption(DICInfo$ExportOptionInfo);
    if(ReadingOptionInfo$GZipped){
      if(is.null(TXTfileName)){
        TXTfileName = paste(.string.cutSubStringFileExtension(DICfileName), ".gz",sep="");
      }
    }else{
      if(is.null(TXTfileName)){
        TXTfileName = paste(.string.cutSubStringFileExtension(DICfileName), ".txt",sep="");
      }
    }
    if("Output filename" %in% DICInfo$SystemInfo$ItemNameInDic){
      DICInfo$SystemInfo["Output filename" == DICInfo$SystemInfo$ItemNameInDic,2] = TXTfileName;
    }else{
      cntRow = .create.emptyDataFrame.string(SectionColNames);
      cntRow[1,1] = "Output filename";
      cntRow[1,2] = TXTfileName;

      cntSectionNameInDIC = attr(DICInfo$SystemInfo, "SectionNameInDIC");
      DICInfo$SystemInfo = rbind(cntRow, DICInfo$SystemInfo);
      attr(DICInfo$SystemInfo, "SectionNameInDIC") = cntSectionNameInDIC; 
    }
    
    UseVarLabelsInDataWhenRead = DICInfo$UseVarLabelsInData;
    if(is.null(UseVarLabelsInDataWhenRead)){
      UseVarLabelsInTxtFile = TRUE;
      DICInfo$VarFormatSecList = list();
    }else{
      if(UseVarLabelsInDataWhenRead){
        if(UseVarLabelsInTxtFile){
          DICInfo$ExportOptionInfo = .assignItemValueToSecInfoVar("quotedlabels","Variable format", DICInfo$ExportOptionInfo); 
          DICInfo$VarFormatSecList = list();
        }
      }else{ #here numeric values in current data object
        if(UseVarLabelsInTxtFile){
        #replace labels for numeric values
          dataInF = .replaceLabOrNum(DICInfo$VarFormatSecList, dataInF, LabReplaced = F)
          DICInfo$ExportOptionInfo = .assignItemValueToSecInfoVar("quotedlabels","Variable format", DICInfo$ExportOptionInfo); 
          DICInfo$VarFormatSecList = list();
        }else{ #do nothing, just output DICInfo$VarFormatSecList and unchanged data object
        
        }
        #UseVarLabelsInTxtFile always true: because data is ready for output, no need of changes
        UseVarLabelsInTxtFile = TRUE;
      }    
    }    
    if(!UseVarLabelsInTxtFile){
      cntItemValue = "numeric";
      DICInfo$ExportOptionInfo = .assignItemValueToSecInfoVar(cntItemValue,"Variable format", DICInfo$ExportOptionInfo); 
    }
  }

  NumVars = length(DICInfo$VarFormatSecList);
  if(!UseVarLabelsInTxtFile && NumVars > 0){
    dataInF = .replaceLabOrNum(DICInfo$VarFormatSecList, dataInF, LabReplaced = T);
  }
  
  attr(dataInF, "DICInfo") = DICInfo;
  return(dataInF);
}


#############################
# exported functions
#############################

#' @export
read.SeerStat = function(DICfileName, TXTfileName=NULL, UseVarLabelsInData=FALSE,ReadHeaderOnly=FALSE,...){
#   input: 
#      DICfileName, TXTfileName 
#   returned: 
#     an object of data.frame storing data in TXTfileName and info in DICfileName
  if(!(length(grep("[.]dic$", DICfileName,ignore.case=TRUE)) == 1)){
    DICfileName = paste(DICfileName, ".dic",sep="");
  }

  return(.read.SeerStatTXT(DICfileName, TXTfileName, UseVarLabelsInData=UseVarLabelsInData,ReadHeaderOnly=ReadHeaderOnly,...));
}


#' @export
write.SeerStat = function(myData, DICfileName, TXTfileName=NULL,UseVarLabelsInTxtFile=TRUE,LabVarsNames=NULL,...){
#   input: 
#      DICfileName, TXTfileName 
#   returned: 
#     a list storing info in DICfileName
  if(!(length(grep("[.]dic$", DICfileName,ignore.case=TRUE)) == 1)){
    DICfileName = paste(DICfileName, ".dic",sep="");
  }

  #prepare DICInfo
  myData = .prepareDICInfoAndData(myData,DICfileName,TXTfileName=TXTfileName,UseVarLabelsInTxtFile=UseVarLabelsInTxtFile,LabVarsNames=LabVarsNames);
  DICInfo = attr(myData, "DICInfo");
  
  .write.SeerStatDIC(DICInfo, DICfileName);
  .write.SeerStatTXT(myData, DICInfo,...);
  return(DICInfo);
}


