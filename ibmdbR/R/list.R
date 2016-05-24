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

################ Constructor ############################
ida.list <- function (type='public',user=NULL) {
  
  idaCheckConnection();
  
  publicTablePrefix <- "R_OBJECTS_PUBLIC";
  privateTablePrefix <- "R_OBJECTS_PRIVATE";
  
  #Get current user
  currentUser <- idaScalarQuery("select trim(current_user) from sysibm.sysdummy1");
  roleName <- idaCheckSharing()
  if(is.null(roleName)&&(type=='public')) {
    warning('Sharing R objects is not supported on this platform, creating private list instead')
    type <- 'private'
  }
  
  #Check if the tables are there already if no other user is given
  if((is.null(user))||(currentUser==user)) {
    if(type=='public') {
      
      tableName <- paste('"',currentUser,'".',publicTablePrefix,sep='')
      tableNameMeta <- paste('"',currentUser,'".',publicTablePrefix,"_META",sep='')
      roleName <- idaCheckSharing()
      if(is.null(roleName)) {
        stop("Sharing R objects is not enabled on this platform.")
      }
      if(!idaExistTable(tableName)) {
          createIdaListTable(tableName,'public',roleName)
      }
      if(!idaExistTable(tableNameMeta)) {
        createIdaListMetaTable(tableNameMeta,'public',roleName)
        createMetaTableContent(tableName,tableNameMeta)
      }
    } else if(type=='private') {
      tableName <- paste('"',currentUser,'".',privateTablePrefix,sep='')
      tableNameMeta <- paste('"',currentUser,'".',privateTablePrefix,"_META",sep='')
      if(!idaExistTable(tableName)) {
        createIdaListTable(tableName,'private')
      }
      if(!idaExistTable(tableNameMeta)) {
        createIdaListMetaTable(tableNameMeta,'private')
        createMetaTableContent(tableName,tableNameMeta)
      }
    } else {
      stop("The only types allowed are public and private")
    }
  } else {
    if(type=='private') {
      stop('You can only access the public folder of other users.')
    } else if(type=='public') {
      tableName <- paste('"',user,'".',publicTablePrefix,sep='')		
    } else {
      stop("The only types allowed are public and private")			
    }
   
    texists <- F
    try({idaQuery("SELECT COUNT(*) FROM ",tableName);texists<-T},silent=T)
    
    if(!texists) {
      stop("The user did not create a table for R object storage yet.")			
    }
  }
  
  return(new(Class="ida.list", tableName=tableName))
}

################ Internal utilities ############################
createIdaListTable <- function(tableName,type,roleName) {

  idaQuery("CREATE TABLE ", tableName, " (OBJID VARCHAR(2000), SID INTEGER, SNIPPET VARCHAR(30000))");
  if(type=='public') {
    idaQuery("GRANT SELECT ON ",tableName," TO ",roleName);
  }
}

createIdaListMetaTable <- function(tableName,type,roleName) {
  idaQuery("CREATE TABLE ", tableName, " (OBJID VARCHAR(2000),SIZE BIGINT, DATE_CREATION VARCHAR(100),CATEGORY VARCHAR(1000))");
  if(type=='public') {
    idaQuery("GRANT SELECT ON ",tableName," TO ",roleName);
  }
}




################ Methods ############################
setMethod("[", signature(x = "ida.list"),
    function(x, i) {
      ida.getObj(i,x@tableName);
    }
)

setMethod("$", signature(x = "ida.list"),
    function(x, name) {
      ida.getObj(name,x@tableName);
    }
)

setMethod("[<-", signature(x = "ida.list"),
    function(x, i,value) {
      ida.storeObj(i,value,x@tableName)
      x;
    }
)

setMethod("$<-", signature(x = "ida.list"),
    function(x, name,value) {
      ida.storeObj(name,value,x@tableName)
      x;
    }
)

setMethod("length", signature(x="ida.list"),
    function(x) { return(idaQuery("select count(distinct OBJID) from ",x@tableName)[[1]][1]) }
)

setMethod("names", signature(x="ida.list"),
    function(x) { return(idaQuery("select distinct OBJID from ",x@tableName)[[1]]) }
)

setMethod("print", signature(x="ida.list"),
    function (x) {
      cat(x@tableName,"\n")
    }
)

setMethod("show", signature(object="ida.list"),
    function (object) {
      cat(object@tableName,"\n")		
    }
)

is.ida.list <- function(x) {
  return(inherits(x, "ida.list"))
}

################ Internal methods for storage and retrieval ############################

ida.storeObj <- function(id, obj, tableName,tableNameMeta=paste(tableName,"_META",sep='')) {
  
  #We need to do this in one transaction
  odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = FALSE)
  
  tryCatch({
        
        #Delete any object with the same key if exists
        ida.delObj(id,tableName);
        
        if(!is.null(obj)) {
          objStr <- rawToChar(serialize(obj, ascii=TRUE,connection=NULL),multiple=FALSE);
          offset <-0;
          maxSnippetSize <- 29000;
          numSnippets <- nchar(objStr)/maxSnippetSize;
          count <- -floor(numSnippets);
          
          while(count<=0){
            objStrSnippet <- substring(objStr,offset,offset+maxSnippetSize-1);
            idaQuery("INSERT INTO ", tableName," VALUES( '", id, "',", count, ",'", objStrSnippet,"')");
            count <- count+1;
            offset<-offset+maxSnippetSize;
          }
          #Get the class and make sure it will fit 
          classStr <- substr(as.character(class(obj)),1,1000);
         
          idaQuery("INSERT INTO ", tableNameMeta," VALUES('", id,"',",max(c(1,as.integer(as.numeric(object.size(obj))/1024.0))),",'", as.character(date()),"','",classStr,"')");
        }
        
        #If all goes well, commit
        odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = TRUE);
      },
      #In cse of error, rollback
      error=function(e){print(e);odbcEndTran(get("p_idaConnection",envir=idaRGlobal), commit = FALSE)},
      #In any case, turn autocommit on again
      finally={odbcSetAutoCommit(get("p_idaConnection",envir=idaRGlobal), autoCommit = TRUE)});
}

ida.delObj <- function(id, tableName,tableNameMeta=paste(tableName,"_META",sep='')) {
  
  try(idaQuery("delete from ",tableName," WHERE OBJID='",id,"'"),silent=T);
  try(idaQuery("delete from ",tableNameMeta," WHERE OBJID='",id,"'"),silent=T);
}

ida.getObj <- function(id, tableName) {
  mdf <- idaQuery("select * from ",tableName," WHERE OBJID='",id,"'");
  
  #Key does not exist
  if(NROW(mdf)==0)
    return(NULL)
  
  obj <- NULL;
  
  df <- mdf[with(mdf, order(SID)), ]
  snippets <- c();
  for(i in 1:nrow(df)) {
    groupid <- df[i,1];
    objStrSnippet<-df[i,3];
    snippetId <- df[i,2];
    snippets <- c(snippets, objStrSnippet);
    if(snippetId == 0) {
      objStr <- paste(snippets,collapse = "",sep="");
      obj<- unserialize(charToRaw(objStr));
    }
  }
  obj
}

createMetaTableContent <- function(snippetTable, metaTable) {
  if(nrow(ida.data.frame(snippetTable))>0) {
    idaQuery("INSERT INTO ",metaTable, " SELECT distinct OBJID, 0, '','' FROM ",snippetTable);
  }
}