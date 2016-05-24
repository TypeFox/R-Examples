import.RInSp = function(filename, col.header=FALSE, row.names = 0, info.cols= 0, subset.column = 0,
                      subset.rows = NA, data.type= "integer", print.messages=TRUE){
    #
    # Work horse function for transforming a resource matrix or dataframe into
    # an object to be used by other functions of the RInSp package
    #
    # Author: Nicola ZACCARELLI
    # E-mail: nicola.zaccarelli@gmail.com
    #
    # Version: 1.0
    # Date: 10/11/2012
    #
    # read data
    if (class(filename) == "character") 
    { if (row.names > 0) {
      datatmp = read.table(filename, header=col.header)
      row.names(datatmp) = datatmp[ , row.names]
    }  
      else datatmp = read.table(filename, header=col.header)
    } 
    else datatmp = filename
    cols = dim(datatmp)[2]
    #
    # some checking before process data
    # check info.cols
    if (info.cols[1] != 0) if (sum(info.cols %in% c(1:dim(datatmp)[2])) != length(info.cols)) stop("Information columns out of range!")
    
    # check subset.column and row names
    if (subset.column[1] != 0) {
      if (sum(subset.column %in% c(1:dim(datatmp)[2])) != length(subset.column)) stop("Wrong subset.column specification. Probably number out of column range.")
      if (row.names %in% subset.column) stop("The row names column must not be part of the subsetting columns set.")}
    # check subset.rows
    if (!is.na(subset.rows[1])) if (length(subset.rows) < 2) stop("Wrong subset.rows specification. There is only one element.") else {
      if (length(grep(subset.rows[1], attributes(datatmp)$names)) == 0) stop("Wrong column's name for rows subsetting.")
      for (i in 2:length(subset.rows)) {
        pos = grep(subset.rows[1], attributes(datatmp)$names)
        if (sum(grep(subset.rows[i], datatmp[ , pos])) == 0) stop("Wrong label for row subsetting.")}
    }
    # check data type
    if (data.type %in% c("integer", "double", "proportion") == FALSE) stop("The specified data type is wrong.")
    
    # subsetting columns and rows
    if (subset.column[1] != 0) column.selection = subset.column else { 
      column.selection = c(1:cols)
      if (info.cols[1] != 0) {
        if (row.names != 0) column.selection = column.selection[-c(row.names, info.cols)]
      } else if (row.names != 0) column.selection = column.selection[-row.names]
    }
    if (!is.na(subset.rows[1])) {
      rows2keep = datatmp[ ,subset.rows[1]]
      for (i in 2:length(subset.rows)){ rows2keep = gsub(subset.rows[i], "XxX", rows2keep) }
      rows2keep = grepl("XxX", rows2keep)
      resources = as.matrix(subset(datatmp, rows2keep, select = column.selection))
      if (info.cols[1] != 0) info = subset(datatmp, rows2keep, select= info.cols) else info = 0
    } else {
      resources = as.matrix(subset(datatmp, select = column.selection))
      if (info.cols[1] != 0) info = subset(datatmp, select= info.cols) else info = 0
    }
    #deleting zero sum rows
    numIndTot = dim(resources)[1]
    numResTot = dim(resources)[2]
    if (info.cols[1] != 0) info = subset(info, apply(resources, 1, sum) > 0)
    resources = subset(resources, apply(resources, 1, sum) > 0) # dropping zero sum diets
    # dropping zero sum resources
    tmp = t(resources)
    tmp = subset(tmp, apply(tmp, 1, sum) > 0) 
    resources = t(tmp)
    numIndEf = dim(resources)[1]
    numResEf = dim(resources)[2]
    row.names = dimnames(resources)[[1]]
    prop = resources / apply(resources, 1, sum)
    col.names = dimnames(resources)[[2]]
    if (data.type == "proportion") ris = list(resources= 0, proportions= prop, data.type= data.type, col.names= col.names, ind.names = row.names, info= info, num.prey= numResEf, num.individuals= numIndEf, num.zero.prey= numResTot - numResEf, num.ind.zero= numIndTot - numIndEf)
    else ris = list(resources= resources, proportions= prop, data.type= data.type, col.names= col.names, ind.names = row.names, info= info, num.prey= numResEf, num.individuals= numIndEf, num.zero.prey= numResTot - numResEf, num.ind.zero= numIndTot - numIndEf)
    if (print.messages == TRUE) {
    if ((numIndTot - numIndEf) != 0) {
      cat("\n Warning! \n")
      cat("\n The total number of sample was", numIndTot, "but", (numIndTot - numIndEf))
      if ((numIndTot - numIndEf) == 1) { cat(" individual was dropped \n")
                                         cat("as not consuming the selected resources \n")} else { 
                                           cat(" individuals were dropped \n")
                                           cat("as not consuming the selected resources \n") }
    }
    if ((numResTot - numResEf) != 0) {
      cat("\n Warning! \n")
      cat("\n The total number of resources was", numResTot, "but", (numResTot - numResEf))
      if ((numResTot - numResEf) == 1) { cat(" resource was dropped \n")
                                         cat("because not present in the selected sample \n")} else {
                                           cat(" resources were dropped \n")
                                           cat("because not present in the selected sample \n") }
    }
	}
    class(ris) = "RInSp"
    return(ris)}

