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
# 

idaMerge <- function(x, y, by = intersect(x@cols, y@cols), by.x = by, by.y = by, all = FALSE, all.x = all, all.y = all, sort = TRUE, suffixes = c("_x", "_y"), table = NULL) {
  
  idaCheckConnection();
  
  #  error check x, y
  if (!is.element(class(x), "ida.data.frame")) {
    stop("idaMerge can only merge two ida.data.frame")
  }
  if (!is.element(class(y), "ida.data.frame")) {
    stop("idaMerge can only merge two ida.data.frame")
  }
  
  #  the initial output columns are the left ON columns
  outputColumnNames = vector()
  outputColumnsForQuery = ""
  orderBy = ""
  remainingColumns1 = x@cols
  remainingColumns2 = y@cols
  for (i in sequence(length(by.x))) {
    if (nchar(outputColumnsForQuery) != 0) {
      outputColumnsForQuery = paste(outputColumnsForQuery, ", ", sep="")
    }
    
    #  notice how we combine the ON columns using nvl, so that the appropriate (non-null) value is shown when doing outer joins
    #  i.e. in R only the single column is outputted whereas in SQL, typically both columns are outputted
    outputColumnsForQuery = paste(outputColumnsForQuery, "nvl(", "x.\"", by.x[i], "\", ", "y.\"", by.y[i], "\") AS \"", by.x[i],"\"", sep="")
    
    outputColumnNames = append(outputColumnNames, by.x[i])
    
    #  remove the left ON column from the columns remaining to be outputted
    remainingColumns1 = remainingColumns1[!toupper(remainingColumns1) == toupper(by.x[i])]
    #  remove the right ON column from the columns remaining to be outputted
    remainingColumns2 = remainingColumns2[!toupper(remainingColumns2) == toupper(by.y[i])]
    
    #  order by the left ON columns
    if (nchar(orderBy) != 0) {
      orderBy = paste(orderBy, ", ", sep="")
    }
    orderBy = paste(orderBy, "x.", by.x[i], sep="")
  }
  
  #  finish creating the order by clause, if any
  if (nchar(orderBy) != 0) {
    orderBy = paste(" ORDER BY ", orderBy, sep="")
  }
  
  
  #  output the remaining left columns
  for (i in sequence(length(remainingColumns1))) {
    if (nchar(outputColumnsForQuery) != 0) {
      outputColumnsForQuery = paste(outputColumnsForQuery, ", ", sep="")
    }
    
    #  if the column name exists append the first suffix (e.g. ".x") to the name
    value = remainingColumns1[i]
    while ((toupper(value) %in% toupper(outputColumnNames)) || (toupper(value) %in% toupper(remainingColumns2))) {
      value = paste(value, suffixes[1], sep="")
    }
    
    outputColumnsForQuery = paste(outputColumnsForQuery, "x.\"", remainingColumns1[i], "\" AS \"", value,"\"", sep="")		
  }
  
  #  output the remaining right columns
  for (i in sequence(length(remainingColumns2))) {
    if (nchar(outputColumnsForQuery) != 0) {
      outputColumnsForQuery = paste(outputColumnsForQuery, ", ", sep="")
    }
    
    #  if the column name exists append the second suffix (e.g. ".y") to the name
    value = remainingColumns2[i]
    while ((toupper(value) %in% toupper(outputColumnNames)) || (toupper(value) %in% toupper(remainingColumns1))) {
      value = paste(value, suffixes[2], sep="")
    }
    
    outputColumnsForQuery = paste(outputColumnsForQuery, "y.\"", remainingColumns2[i], "\" AS \"", value, "\"",sep="")		
  }
  
  #  generate ON clause from the by parameters
  onClause = ""
  for (i in sequence(length(by.x))) {
    if (nchar(onClause) != 0) {
      onClause = paste(onClause, " AND ", sep="")
    }
    onClause = paste(onClause, "x.\"", by.x[i], "\" = ", "y.\"", by.y[i],"\"", sep="")
  }
  #  if no ON clause, then this will be a cartesian product
  if (nchar(onClause) == 0) {
    joinType = " CROSS JOIN "
  }
  #  else there is an ON clause
  else {
    onClause = paste(" ON ", onClause, sep="")
    
    
    if (!all && !all.x && !all.y) {
      joinType = " INNER JOIN "
    }
    else if (all.x && all.y){
      joinType = " FULL OUTER JOIN "
    }
    else if (all.x) {
      joinType = " LEFT OUTER JOIN "
    }
    else {
      joinType = " RIGHT OUTER JOIN "
    }
  }
  
  #  if no table name provided, then generate a unique one
  if (is.null(table))
    table = idaGetValidTableName()
  
  # If an input table contains a where clause, we create a view first, before
  if(is.null(x@where)&&(nchar(x@where)>0)) {
    xView <- idadf.from(x);
  } else {
    xView <- paste("(",idadf.query(x),")",sep='');
  }
  
  if(is.null(y@where)&&(nchar(y@where)>0)) {
    yView <- idadf.from(y);
  } else { 
    yView <- paste("(",idadf.query(y),")",sep='');
  }
  
  query = paste("create view ", table, " as select ", outputColumnsForQuery, " from ", xView, " x", joinType, yView, " y", onClause, sep = "")
  idaQuery(query)
  
  output = ida.data.frame(table)
  
  return(output)
}

