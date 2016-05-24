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


idaInit <- function(con,jobDescription=NULL) {
  #Check if the connection is open?
  conOpen <- FALSE;
  try({idadf(con,'select trim(current_schema) from sysibm.sysdummy1');conOpen<-TRUE},silent=TRUE);
  if(!conOpen) {
    stop("con is not an open connection, please use idaConnect() to create an open connection to the data base.");
  }
  
  #Put the connection object into a new environment
  assign("idaRGlobal", new.env(parent = baseenv()), envir=baseenv())
  assign("p_idaConnection", con, envir = idaRGlobal)
  
  #Check Oracle compatibility
  regVars <- idaQuery("SELECT reg_var_name, reg_var_value, level FROM table(REG_LIST_VARIABLES()) as reg")
  compVec <- regVars[regVars$REG_VAR_NAME=='DB2_COMPATIBILITY_VECTOR',2] 
  
  if(length(compVec)) {
    if(compVec=='ORA') {
      assign("p_db2_comp", "ORA", envir = idaRGlobal) 
    }
  }
  
  script <- jobDescription
  #Set the application parameters
  if(is.null(script)) {
    #Need to find out how the current script is called if possible
    script <- 'interactive' 
    scriptFiles <- lapply(sys.frames(), function(x) x$ofile)
    for(i in 1:length(scriptFiles)) {
      if(!is.null(scriptFiles[[i]])) {
        script <- scriptFiles[[i]];
        break;
      }
    }
  }
    
  #Set the current connection parameters for job monitoring
  try({idaQuery(paste("CALL WLM_SET_CLIENT_INFO(NULL,NULL,'ibmdbR','",script,"',NULL)",sep=''))},silent=TRUE)
   
  #Check what functions are available in the database
  c1 <- idaCheckProcedure("KMEANS","idaKMeans",TRUE)
  c2 <- idaCheckProcedure("NAIVEBAYES","idaNaiveBayes",TRUE)
  c3 <- idaCheckProcedure("ASSOCRULES","idaArule",TRUE)
  c4 <- idaCheckProcedure("LINEAR_REGRESSION","idaLm",TRUE)
  #c5 <- idaCheckProcedure("SEQRULES","idaSeqRules",TRUE)
  c5 <- T
  
  
  if(!(c1&&c2&&c3&&c4&&c5)) {
    message("Note that not all backend databases provide push-down capabilities for all analytical functions.")
  }
  
}

idaIsOracleMode <- function() {
  return(exists("p_db2_comp",envir=idaRGlobal)&&(get("p_db2_comp",envir=idaRGlobal)=='ORA'))
}

idaCheckConnection <- function () {
  # make sure that the connection exists
  if ((!exists("p_idaConnection",envir=idaRGlobal)) || is.null(get("p_idaConnection",envir=idaRGlobal)))
    stop("The connection is not set, please use idaInit(con), where con is an open connection.", call.=FALSE)
}

idaCheckProcedure <- function(procName, algName, verbose=FALSE) {
  
  catQuery <- paste("SELECT COUNT(*) FROM SYSCAT.ROUTINES WHERE ROUTINENAME = '",procName,"' AND ROUTINEMODULENAME = 'IDAX'",sep='') 
  available <- as.numeric(idaScalarQuery(catQuery))>0;
  
  if(verbose&&!available) {
    message(paste("Function ",algName, " ",ifelse(available,"","not "),"available for this connection.",sep=''))
  } 
  
  invisible(available)
}

idaCheckRole <- function(roleName) {
  #catQuery <- paste("SELECT COUNT(*) FROM SYSCAT.ROLES WHERE ROLENAME = '",roleName,"'",sep='') 
  #available <- idaScalarQuery(catQuery)>0;
  
  #if(!available)
  #  return(FALSE)
  
  catQuery <- paste("SELECT COUNT(*) FROM SYSCAT.ROLEAUTH WHERE ROLENAME = '",roleName,"' AND GRANTEE=CURRENT USER",sep='') 
  granted <- as.numeric(idaScalarQuery(catQuery))>0;
  
  return(granted)
}

idaCheckSharing <- function() {
  
  roledashDB <- "DASHDB_ENTERPRISE_USER"
  roleDB2 <- "R_USERS_PUBLIC"

  #First try dashDB
  if(idaCheckRole(roledashDB)) {
    return(roledashDB)
  } else if(idaCheckRole(roleDB2)) {
    return(roleDB2)
  } else {
    return(NULL)
  }
  
  
}

