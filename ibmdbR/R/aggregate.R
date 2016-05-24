# 
# Copyright (c) 2010, 2014, IBM Corp. All rights reserved. 
# 		
# This program is free software: you can redistribute it and/or modify 
# it under the terms of the GNU General Public License as published by 
# the Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version. 
#
# This program is distributed in the hope that it will be useful, 
# but WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU General Public License for more details. 
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>. 
#
idaAggregate <- function(data,INDICES,FUN,output.name=NULL,clear.existing=F) {
  
  if (!is.ida.data.frame(data))
    stop('data has to be a ida.data.frame')
  
  aggMap <- list()
  aggMap["mean"] <- "AVG"
  aggMap["min"] <- "MIN"
  aggMap["max"] <- "MAX"
  aggMap["var"] <- "VARIANCE"
  aggMap["sum"] <- "SUM"
  aggMap["sd"] <- "STDDEV"
  aggMap["length"] <- "COUNT"

  if((!is.list(INDICES))&&(!inherits(INDICES,"ida.col.def")))
    stop("INDICES must either be a list or a single column definition.")
  
  #First try whether the function is specified as string
  if(!is.character(FUN)) {
    #Now try whether the function symbol name is a function we know
    FUN <-  as.character(substitute(FUN))
  }
  print(FUN)  
  sqlFun <- aggMap[[FUN]]
  if(is.null(sqlFun)) {
    stop("Function not supported, currently only mean, min, max, var, sum, sd and length are supported.")
  }
  
  indexTerms <- c()
  indexCols <- c()
  if(is.list(INDICES)) {
    
    #Add column names if needed
    if(is.null(names(INDICES))) {
      names(INDICES) <- as.character(1:length(INDICES))
    }
    
    for(i in 1:length(names(INDICES))) {
      #Check column type
      name <- names(INDICES)[i]
      print(name)
      q <- INDICES[[name]];
      print(q)
      if(!inherits(q,"ida.col.def")) {
        stop("INDICES need to be columns of an ida.data.frame, e.g. idf$X")
      }
      indexTerms <- c(indexTerms,q@term)
      indexCols <- c(indexCols,name)
      #Add to data column and coldef
    }
    groupByPart <- paste(indexTerms,collapse=",",sep="")
    selectPart <- paste(indexTerms,' AS "',indexCols,'"',collapse=",",sep="")
  } else {
    groupByPart <- INDICES@term
    selectPart <- paste(INDICES@term,' AS "',c(substitute(INDICES)),'"',sep="")
  }
  
  query <- paste("SELECT ", selectPart,",",paste(sqlFun,'("',data@cols,'") AS "',sqlFun,'(',data@cols,')"',collapse=",",sep='')," FROM ", idadf.from(data), ifelse((nchar(data@where)>0), paste(" WHERE ", data@where), ""), "GROUP BY ",groupByPart);
  print(query)	
  
  if(!is.null(output.name)) {
    if(idaExistTable(output.name)) {
      if(clear.existing) {
        idaDeleteViewOrTable(output.name)
      } else {
        stop("View already exists, use the clear.existing option to delete.")
      }
    }
    
    idaQuery("CREATE VIEW ",output.name," AS ",query);
    return(ida.data.frame(output.name))
  } else {
    return(idaQuery(query))
  }
}
