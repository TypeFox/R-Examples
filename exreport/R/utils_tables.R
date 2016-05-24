

.orderCascade <- function(df, decreasing = NULL){
  # Order a data.frame in cascade, by columns order (for each column the vector decreasing determines the order)
  # Non-numeric columns are ignored.
  # Returns the correct index order for the data.frame
  
  # Only numeric columns are of interest.
  colOfInterest <- sapply(df,is.numeric)
  numericDf <- df[colOfInterest]
  numericDecreasing <- decreasing[colOfInterest]
  # Matrix to make the opposite if decreasing equals TRUE. (each column have
  # one unique value. +1 if the order is increasing and -1 if it is decreasing). 
  orderingMatrix <- matrix(rep(numericDecreasing*-2+1,nrow(numericDf)),nrow = nrow(numericDf),byrow = TRUE)
  reorderedMatrix <- numericDf*orderingMatrix
  colnames(reorderedMatrix) <- NULL
  # Get the ordered indices
  ordering <- do.call(function(...) order(...,na.last = FALSE), lapply(reorderedMatrix,FUN = identity))
  ordering
}


.duplicates <- function(dat){
  #Search for duplicated rows in a data.frame and return the index for that occurrences. NAs for the first occurrences.
  #
  # Args:
  #   dat:           The data.frame
  #
  # Returns:
  #   a numeric array (the index of the original rows if duplicated or NA).
  
  # PARAMETER VALIDATION:
  if (!is.data.frame(dat))
    stop("error, you must provide a data.frame object as 'dat' parameter")
  # Get the indices of the ordered rows
  unamedDat <- dat
  colnames(unamedDat) <- NULL
  s <- do.call("order",lapply(unamedDat,FUN=identity))
  # As we use the ordered dat, for each set of duplicated rows, the first occurrence is TRUE and the laters are FALSE
  nonDup <- !duplicated(dat[s, ,drop=FALSE])
  # Get the index of the non duplicated (or original) elements (elements equal to TRUE)
  origInd <- s[nonDup]
  # cumsum applies the cumulative sum for the nonDup (casting logical to numeric, so 0 and 1).
  # This way, for each set of duplicated rows, all the elements will have the same number (only adds 1 the first element because it is TRUE)
  firstOcc <- origInd[cumsum(nonDup)]
  # We set to NA the elements that are non duplicated (or original).
  firstOcc[nonDup] <- NA
  # The rest of the elements will have the index of his original row (as we unorder the firstOcc array, we get the original order --> the same as dat)
  firstOcc[order(s)]
} 

.formatDataFrame <- function(tables, formats, src) {
  # df: the data.frame
  # src: latex/html
  # formats: vector of formats for the data.frame columns (in the sprintf format)
  # thresholds: vector of fix/max/min
  # values: vector of values for the function
  if (src == "html")
    res <- .formatDataFrameHTML(tables, formats)
  else if (src == "latex")
    res <- .formatDataFrameLatex(tables, formats)
  else
    stop("Invalid src")
  res
}

.formatDataFrameHTML <- function(tables, formats, charForNAs = "-") {
  cnames <- colnames(tables[[1]])
  # We apply the format process to all tables in tables, and then we put this information together to build the final HTML table.
  strFormattedTables <- list()
  # Function to translate a format data.frame from "b" to "<strong>s/strong>"
  # Now only b (bold) is allowed. If we want to add more (i for italic i.e.), is just the same:
  # extract the logical array, remove the letter, and apply the ifelse.
  funTranslateFormat <- function(format){
    apply(format,2,FUN=function(v){
      isStrong <- grepl("b",v)
      aux <- gsub("b","",v)
      ifelse(isStrong,paste0("<strong>",aux,"</strong>"),aux)
    })
  }
  # Function to mix two data.frames (table and format) into a matrix (of strings) according to sprintf syntax.
  funApplyFormat <- function(table, format){
    strTable <- matrix("",nrow = nrow(table),ncol=ncol(table))
    for(i in 1:ncol(table)){
      strTable[,i] <- ifelse(is.na(table[,i]),charForNAs,sprintf(format[,i],table[,i]))
    }
    strTable
  }
  
  
  # Function to translate from a matrix of strings to a single string in tabular HTML syntax
  funMatrix2Str <- function(m){
    apply(m, 1, FUN=function(x){
      sprintf("<tr>\n%s%s\n</tr>\n",sprintf("<td style=\"text-align:left\">%s</td>",x[1]), paste0(sprintf("<td style=\"text-align:right\">%s</td>",x[-1]),collapse = ""))
    })
  }
  
  # For each table, we apply the funTranslateFormat function
  auxFormats <- lapply(formats, FUN = funTranslateFormat)
  # For each pair table-format, we apply the funApplyFormat function
  auxMatrix <- mapply(FUN = funApplyFormat, tables, auxFormats, SIMPLIFY = FALSE)
  
  # Now we have final HTML tables for individual outputs. We have to combine them into a single one.
  # We first calculate the desired indices for the matrix when combine all of them into a single one.
  colIndices <- rep(2:ncol(tables[[1]]),each=length(tables)) + (ncol(tables[[1]])-1)*(seq(1,length(tables))-1)
  allMatrix <- auxMatrix[[1]]
  if(length(auxMatrix)>1){
    for(i in 2:length(auxMatrix))
      allMatrix <- cbind(allMatrix,auxMatrix[[i]][,-1,drop=FALSE])
  }
  allMatrix <- allMatrix[,c(1,colIndices)]
  
  # Now we translate from matrix to tabular HTML syntax
  strTables <- funMatrix2Str(allMatrix)
  
  # Finally, we build the entire latex table
  res <- "<table class=\"table table-hover\">\n<tr>\n"
  if(length(tables)==1)
     res <- paste0(res, "<tr>\n<th style=\"text-align:left\">",cnames[1],"</th>", paste0(sprintf("<th style=\"text-align:right\">%s</th>",cnames[-1]),collapse = ""), "</tr>\n")
  else{
    res <- paste0(res, "<th></th>", paste0(sprintf("<th style=\"text-align:center\" colspan=\"%d\">%s</th>",length(tables),cnames[-1]),collapse = ""), "</tr>\n")
    #Ahora los outputs
    res <- paste0(res, "<tr>\n<th style=\"text-align:left\">",cnames[1],"</th>",paste0(sprintf("<th style=\"text-align:right\">%s</th>",paste0(rep(names(tables), length(cnames)-1)),collapse = ""),collapse = ""),"</tr>\n")
  }
  
  res <- paste0(res,paste0(strTables, collapse = ""))
  
  res <- paste(res,"</table>",sep="")
  res
}

.formatDataFrameLatex <- function(tables, formats, charForNAs = "-") {
  cnames <- colnames(tables[[1]])
  # We apply the format process to all tables in tables, and then we put this information together to build the final Latex table.
  strFormattedTables <- list()
  # Function to translate a format data.frame from "b" to "\\bf
  # Now only b (bold) is allowed. If we want to add more (i for italic i.e.), is just the same:
  # extract the logical array, remove the letter, and apply the ifelse.
  funTranslateFormat <- function(format){
    apply(format,2,FUN=function(v){
      isBf <- grepl("b",v)
      aux <- gsub("b","",v)
      ifelse(isBf,paste0("\\bf ",aux),aux)
    })
  }
  # Function to mix two data.frames (table and format) into a matrix (of strings) according to sprintf syntax.
  funApplyFormat <- function(table, format){
    strTable <- matrix("",nrow = nrow(table),ncol=ncol(table))
    for(i in 1:ncol(table)){
      strTable[,i] <- ifelse(is.na(table[,i]),charForNAs,sprintf(format[,i],table[,i]))
    }
    strTable
  }
  
  
  # Function to translate from a matrix of strings to a single string in tabular latex syntax
  funMatrix2Str <- function(m){
    apply(m, 1, FUN=function(x){
      paste0(paste(x, collapse=" & "),"\\\\\n")
    })
  }
  
  # For each table, we apply the funTranslateFormat function
  auxFormats <- lapply(formats, FUN = funTranslateFormat)
  # For each pair table-format, we apply the funApplyFormat function
  auxMatrix <- mapply(FUN = funApplyFormat, tables, auxFormats, SIMPLIFY = FALSE)
  
  # Now we have final latex tables for individual outputs. We have to combine them into a single one.
  # We first calculate the desired indices for the matrix when combine all of them into a single one.
  colIndices <- rep(2:ncol(tables[[1]]),each=length(tables)) + (ncol(tables[[1]])-1)*(seq(1,length(tables))-1)
  allMatrix <- auxMatrix[[1]]
  if(length(auxMatrix)>1){
    for(i in 2:length(auxMatrix))
      allMatrix <- cbind(allMatrix,auxMatrix[[i]][,-1,drop=FALSE])
  }
  allMatrix <- allMatrix[,c(1,colIndices)]
  
  # Now we translate from matrix to tabular latex syntax
  strTables <- funMatrix2Str(allMatrix)


  # Finally, we build the entire latex table
  res <- "\\begin{tabular}{l"
  res <- paste0(res, paste0(rep("r",ncol(allMatrix)-1), collapse=""), "}\n")
  if(length(tables)==1)
    res <- paste0(res, "\\hline\n", paste(cnames, collapse = " & "), "\\\\\n\\hline\n")
  else{
    res <- paste0(res, "\\hline\n")
    res <- paste0(res, " & ", paste0("\\multicolumn{",length(tables),"}{c}{",cnames[-1],"}", collapse = " & "), "\\\\\n")
    #Ahora los outputs
    res <- paste0(res, cnames[1], " & ",paste0(rep(names(tables), length(cnames)-1),collapse = " & "), "\\\\\n\\hline\n")
  }
  
  res <- paste0(res,paste0(strTables, collapse = ""))
  
  res <- paste0(res,"\\hline\n\\end{tabular}\n")
  res
}

.mlTableTranspose <- function(exTabular, firstColName){
  #Transpose both table and format data.frames. Column names will be as the first column, and the
  newTables <- list()
  tables <- exTabular$tables
  for(i in 1:length(tables)){
    #first column will be as columnames
    newTable <- tables[[i]]
    #We remove the first column
    fc <- newTable[,1]
    newTable <- newTable[,-1]
    newTable <- t(newTable)
    #Now set the rownames as the first column
    newTable <- data.frame(rownames(newTable),newTable,stringsAsFactors = FALSE)
    rownames(newTable)<-NULL
    #Now set as colnames the new firstColName and fc
    colnames(newTable) <- c(firstColName,as.character(fc))
    newTables[[names(tables)[i]]] <- newTable
  }
  #The same operations for the data.frame format
  newFormats <- list()
  formats <-  exTabular$formats
  for(i in 1:length(formats)){
    newFormat <- formats[[i]]
    #We remove the first column
    newFormat <- newFormat[,-1]
    newFormat <- t(newFormat)
    #Now set as the first column all %s
    newFormat <- data.frame(rep("%s",nrow(newFormat)),newFormat,stringsAsFactors = FALSE)
    rownames(newFormat)<-NULL
    colnames(newFormat) <- colnames(newTables[[i]])
    newFormats[[names(formats)[i]]] <- newFormat
  }
  
  tab <- .exTabular(tables = newTables,
                    formats=newFormats,
                    tableSplit=exTabular$tableSplit,
                    tableType= exTabular$tableType,
                    title = exTabular$title,
                    tags = exTabular$tags)
  tab
}