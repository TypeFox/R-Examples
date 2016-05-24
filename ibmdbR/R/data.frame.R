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

relop <- data.frame(op=c("<",">", "<=", ">=", "!=", "=="),
    sql=c("","","","","<>", "="))

logop <- data.frame(op=c("&", "|"), sql=c("AND", "OR"))

################ Constructor ############################

ida.data.frame <- function (table) {
  
  if(!idaExistTable(table)) {
    stop("Table ", table, " does not exist");
  }
  
  cols <- idaListTableColumns(table)
  
  for (i in 1:length(cols)) 
    if (nchar(cols[i]) > 2 && substr(cols[i], 1, 1) == '"' && substr(cols[i], nchar(cols[i]), nchar(cols[i])) == '"') {
      
      #remove surrounding \"
      cols[i] <- substr(cols[i],2,nchar(cols[i])-1)
    }
  
  return(new(Class="ida.data.frame", table=table  , where="", cols=cols,colDefs=list()));
}

idadf.defined.columns <- function(idadf) {
  
  cols <- idadf@cols
  defCols <- c()
  for(i in 1:length(cols)) {
    if(!is.null(idadf@colDefs[[cols[i]]])) {
      defCols <- c(defCols,cols[i])
    }
  }
  defCols
}

################ Query and internal utilities ############################

idadf.query <- function (idadf, max.rows=NULL, order.by=NULL,includeWhere=T) {
  
  if (!is.ida.data.frame(idadf))
    stop(paste("cannot query object of class: ", class(idadf)))
  
  cols <- idadf@cols
  colTerms <- c();
  for(i in 1:length(cols)) {
    if(is.null(idadf@colDefs[[cols[i]]])) {
      colTerms <- c(colTerms,paste('"', cols[i],'"', sep=''))
    } else {
      colTerms <- c(colTerms,paste('(',idadf@colDefs[[cols[i]]],') AS "', cols[i],'"', sep=''))
    }
  }
  
  paste("SELECT ", paste( colTerms, collapse=","), " FROM ",
      idadf@table,
      ifelse((nchar(idadf@where)>0)&&includeWhere, paste(" WHERE ", idadf@where), ""),
      ifelse(is.null(order.by), "", paste(" ORDER BY ", order.by, sep="")),
      ifelse(is.null(max.rows), "", paste(" FETCH FIRST ", max.rows, " ROWS ONLY")),
      sep='')
}

idadf.from <- function (x) {
  if (!is.ida.data.frame(x))
    stop('x has to be a ida.data.frame')
  
  if(length(idadf.defined.columns(x))>0) {
    x@cols <- union(x@cols, idaListTableColumns(x@table))
    return(paste("(",idadf.query(x,includeWhere = F),")",sep=''))
  } else {
    return(x@table);
  }
}

is.ida.data.frame <-
    function(x) {
  return(inherits(x, "ida.data.frame"))
}

################ Generic methods for upload and download of data ############################

# ---------------------------------------------------------------------

setMethod("as.data.frame", signature(x="ida.data.frame"), 
    function (x, max.rows=NULL, ...) {
      #TODO(opt): BIGINT automatic conversion
      ans <- idaQuery(idadf.query(x, max.rows=max.rows))	
      return(ans)
    }
)

# ---------------------------------------------------------------------
as.ida.data.frame <- function (x, table = NULL, clear.existing=FALSE, case.sensitive=TRUE) {
  
      if(!is.null(table)) {
        
        if(toupper(table) != table)
          warning("Mixed case table names are not supported currently.")
        
        tabSchema <- parseTableName(table,FALSE,TRUE);
        table <- toupper(paste(tabSchema$schema,".",tabSchema$table,sep=''));
      }
      
      if(!is.null(table)&&idaExistTable(table)) {
        if(clear.existing) {
          idaDeleteTable(table);
        } else {
          stop("Table already exists, choose different name or use clear.existing option.")	
        }
      }	
      
      if(is.null(table)) {
        table <- idaGetValidTableName();
      } 
      
      if(!case.sensitive) {
        names(x) <- toupper(names(x));
      }
      
      sqlSave(get("p_idaConnection",envir=idaRGlobal), dat=x, tablename = table, rownames=F,fast=ifelse(idaIsOracleMode(),F,T));
      return(ida.data.frame(as.character(table)));		
    }


################ Sub data frames ############################

setMethod("[", signature(x = "ida.data.frame"), 
    function (x, i=NULL, j=NULL, ..., drop=NA)
    {
      c <- c()
      # check arguments - columns
      if (try(!is.null(j),silent=TRUE) == TRUE) {
        if (is.numeric(j)){
          colCount <- length(x@cols)
          for (n in j) {
            pos <- which(j==n)
            if (n > colCount) {
              stop(paste('The ', pos, '-th column index ', n, 
                         ' in the argument is out of boundary, the maximal column index is ', 
                         colCount, sep=''))
            }
            if (length(pos)>1) {
              stop(paste('The column index ', n, ' appears more than once in the argument at position ', 
                         paste(pos, collapse=', '), '.', sep=''))
            } 
          }
          c <- c(c,as.integer(j))
        }
        else if (is.logical(j)){
          if (length(j) > length(x@cols))
            length(j) <- length(x@cols)
          c <- c(c,which(j))
        }
        else if (!is.integer(j)){
          if (is.character(j)){
            for (n in j){
              if (is.element(n, x@cols))
                c <- c(c, which(names(x)==n))
              else if (is.element(tolower(n), x@cols))
                stop(paste("No column named ", n, " in the table. Column names are case-sensitive. Did you mean ", tolower(n), "?"))
              else if (is.element(toupper(n), x@cols))
                stop(paste("No column named ", j, " in the table. Column names are case-sensitive. Did you mean ", toupper(n), "?"))
              else 
                stop(paste("No column named ", n, " in the table. Column names are case-sensitive.  Candidates are: ", paste(idaListTableColumns(x@table), collapse=', '), "."))
              pos <- which(j==n)
              if (length(pos)>1) {
                stop(paste("The name '", n, "' appears more than once in the column argument at position ", 
                           paste(pos, collapse=', '), '.', sep=''))
              }
            }
          }
          else
            stop("columns argument must be integer or character")
        }
      }
      # check arguments i (row)
      
      if(!missing(i)) {
        if (tryCatch(!is.null(i),error = function(e) {print(e);print("Sub set selection could not be created, the left-hand side of the expression must be a column reference, the right-hand side must be a value or a column reference in the same table.")}) == TRUE) {
          if (is.numeric(i))
            stop("row numbering is not allowed")
          else if (class(i)=="ida.col.def") {
            if((i@table@table != x@table)||(i@table@where!=x@where))
              stop("Cannot apply condition to columns not in the base table.")
            
            if(i@type!='logical')
              stop("Column expression must resolve into a boolean value for row selection.")            
            newRowSelection <- i@term;
          } 
          else if (class(i) == "ida.data.frame.rows") 
            newRowSelection <- i@where
          else 
            stop("row object does not specify a subset")
          
          if (is.null(x@where) || !nchar(x@where))
            x@where <- newRowSelection
          else
            x@where <- paste("(", x@where, ") AND (", newRowSelection, ")", sep="")
        }
      }
      # compute the right subset of columns
      if (!is.null(x@cols) && !is.null(c))
        x@cols <- x@cols[c]
      # i variable has to be of a special class "ida.data.frame.rows"
      return(x)
    }
)

# ---------------------------------------------------------------------

setMethod("$", signature(x = "ida.data.frame"),
    function(x, name) {
      
      if(!(name %in% x@cols))
        stop("Column not found in ida.data.frame.")
      
      if(is.null(x@colDefs[[name]])) {
        return(new(Class="ida.col.def",term=paste('"',name,'"',sep=''),table=x,type="column",aggType="none"));
      } else {
        return(new(Class="ida.col.def",term=x@colDefs[[name]],table=x,type="expr",aggType="none"));
      }
    }
)

setMethod("$<-", signature(x = "ida.data.frame"),
    function(x, name,value) {
      if(is.null(value)) {
        #remove col def
        if(!is.null(x@colDefs[[name]])) {
          x@colDefs[[name]] <- NULL;
        } 
        x@cols <- setdiff(x@cols,name)
      } else {
        
        if(!inherits(value,"ida.col.def"))
          stop("Column definition is not valid for a ida.data.frame, please refer to the documentation of ida.data.frame for details on usage.")
        
        if(((value@table)@table != x@table)||((value@table)@where!=x@where))
          stop("Column definitions are only allowed on the same underlying table. Please use idaMerge first to join the tables.");
        
        if(value@aggType!='none') {
          stop("Cannot add column that contains aggregation term to ida.data.frame.")
        }
        
        x@colDefs[[name]]<-value@term;
        
        if(!(name %in% x@cols)) {
          x@cols<-c(x@cols,name);
        }
      }
      
      return(x);
    }
)

idaCreateWherePart <- function (obj, value, operator) {
  cols <- obj@cols
  
  return(new(Class="ida.data.frame.rows",
          where = paste('"',cols, "\" ", operator, " ", value,  sep="", collapse=" AND ")))
}

invisible(apply(relop, 1,
        function(x) setMethod(x[1], signature(e1="ida.data.frame"),
              function (e1, e2) {
                
                #We allow for column references on the right hand side but only if a single column is referenced and if both data frames point
                #to the same table in the data base
                if(inherits(e2,'ida.data.frame')) {
                  if((length(e2@cols)==1)&&(e1@table==e2@table)) {
                    value <- paste('"',e2@cols[1],'"',sep='');
                  } else {
                    value <- e2;
                  }
                } else {
                  value <- paste("'",e2,"'",sep='');
                }
                
                return(idaCreateWherePart(e1, value, ifelse(nchar(x[2]), x[2], x[1])))
              }
          )
    ))


invisible(apply(logop, 1,
        function(x) setMethod(x[1],
              signature(e1="ida.data.frame.rows",	e2="ida.data.frame.rows"),
              function(e1, e2)
                return (new(Class="ida.data.frame.rows",
                        where=paste("(", e1@where, ") ", x[2], " (", e2@where, ")")))
          )
    ))


################ Basic ida.data.frame functions and methods ############################

setMethod("print", signature(x="ida.data.frame"),
    function (x) {
      cat(idadf.query(x),"\n")
    }
)

setMethod("show", signature(object="ida.data.frame"),
    function (object) {
      cat(idadf.query(object),"\n")
    }
)

idaCreateView <- function (x) {
  
  idaCheckConnection();
  
  if (!is.ida.data.frame(x))
    stop("idaCreateView is valid only for ida.data.frame objects")
  
  name <- idaGetValidTableName("idar_view_")
  idaQuery("CREATE VIEW ", name, " AS ", idadf.query(x))
  return(name)
}
