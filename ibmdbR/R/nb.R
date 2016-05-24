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

idaNaiveBayes <- function(form, data, id="id", modelname=NULL) {
  
  if(!idaCheckProcedure("NAIVEBAYES","idaNaiveBayes",F)) {
    stop("Function not available.")
  }
  
  ntab <- idaParseRFormula(form, data)
  varY  <- paste('\"',ntab$response,'\"',sep="")
  
  model <- modelname
  if (is.null(model)) {
    model <- idaGetValidModelName('NB_')
  } else {
    if(grepl(" ",model)) {
      stop("Space in model name not allowed.")
    }
    
    xx <- parseTableName(modelname);
    model <- paste('"',xx$schema,'"."',xx$table,'"',sep=''); 	
  }
  
  
  colu <- data@cols
  if (!(id %in% colu))
    stop(simpleError(paste("Id variable is not available in ida.data.frame:", id)))
  
  dataTmp <- data[,which(data@cols %in% c(ntab$cols, id,ntab$response))]
  tmpView <- idaCreateView(dataTmp)
  
  id  <- paste('"',id,'"',sep='');
  
  tryCatch({
        callSP("IDAX.NAIVEBAYES ", model=model, intable=tmpView, id=id, target=varY)
      }, error = function(e, tmpView) {
        # in case of error, drop view and let user know, what happend
        stop(e)
      }, finally = {
        # drop view
        idaDropView(tmpView)
      }
  )
 
  result <- idaRetrieveNBModel(model)
  
  return(result)
}

#------------------------------------------------------------------------------

idaRetrieveNBModel <- function(modelName) {
  
  xx <- parseTableName(modelName);
  model <- xx$table
  modelSchema <- xx$schema
  
  
  # results are converted to naiveBayes object
  model2 <- paste('SELECT * FROM "',modelSchema,'"."',model,'_MODEL"',sep="")
  nbOut <- idaQuery(model2);
  
  nbOut$CLASSVALCOUNT <- as.numeric(nbOut$CLASSVALCOUNT)
  nbOut$ATTRCLASSCOUNT <- as.numeric(nbOut$ATTRCLASSCOUNT)
  nbOut$CLASSCOUNT <- as.numeric(nbOut$CLASSCOUNT)
  nbOut$TOTALCOUNT <- as.numeric(nbOut$TOTALCOUNT)
  
  class <- "CLASS"
  val <- "VAL"
  attribute <- "ATTRIBUTE"
  classvalcount <- "CLASSVALCOUNT"
  attrclasscount <- "ATTRCLASSCOUNT"
  colname <- "COLNAME"
  totalcount <- "TOTALCOUNT"
  
  classes <- unique(nbOut[[class]])	
  tables <- list()
  xVars <- unique(nbOut[[attribute]])
  for (i in 1:length(xVars)) {
    xVar <- xVars[i]
    tmp <- nbOut[nbOut[[attribute]]==xVar, ]
    values <- unique(tmp[[val]])
    tab <- matrix(0, length(classes), length(values))
    names.list <- list()
    names.list[["Y"]] <- classes
    names.list[[xVar]] <- values
    dimnames(tab) <- names.list
    for (j in 1:nrow(tmp)) {
      tab[tmp[j,class], tmp[j,val]] <- as.numeric(tmp[j,classvalcount])/tmp[j,"CLASSCOUNT"]
    }
    tables[[i]] <- tab
    names(tables)[i] <- xVar
  }
  splits.rcv <- NULL
  # check if there are automatically made splits for continuous values
  disc.name <- paste('"',modelSchema,'"."',model,'_DISC"',sep="")
  if (idaExistTable(disc.name)) {
    splits.rcv <- idaQuery(paste('SELECT * FROM "',modelSchema,'"."',model,'_DISC"',sep=""))
    classes <- unique(splits.rcv[[colname]])
    splits <- list()
    for (i in 1:length(classes)) {
      splits[[classes[i]]] <- splits.rcv[splits.rcv[[colname]] == classes[i], "BREAK"]
    }
  } else {
    splits <- list()
  }
  
  # apriori property
  apriori <- tapply(nbOut[,6], nbOut[,3], mean, na.rm=T)[unique(nbOut[[class]])]
  
  # prepare object of the class idaNaiveBayes
  new.nb <- list(apriori = apriori, tables = tables, levels = classes, splits = splits, call = "n/a", model = model, modelTable = model2)
  
  class(new.nb) <- c("idaNaiveBayes", "naiveBayes")
  
  return(new.nb)
}

print.idaNaiveBayes <- function (x, ...) {
  cat("Naive Bayes model:\n")
  cat("Apriori probabilities:\n")
  
  n <- sum(x$apriori)
  aprioriProb <- x$apriori/n
  
  print(aprioriProb, ...)
  
  cat("\nConditional probabilities:\n")
  print(x$tables, ...)
  invisible(x)
}


predict.idaNaiveBayes <- function(object, newdata, id,...) {
  
  newData <- newdata
  
  outtable <- idaGetValidTableName(paste("PREDICT_",sep=""))
  
  colu = newData@cols
  if (!(id %in% colu))
    stop(simpleError(paste("Id variable is not available in ida.data.frame:", id)))
  
  id  <- paste('\"',id,'\"',sep="")
  tmpView <- idaCreateView(newData)
  
  tryCatch({	
          callSP("IDAX.PREDICT_NAIVEBAYES ",
              model=object$model,
              intable=tmpView,
              id=id,
              outtable=outtable,
              ...);
          
      }, error = function(e) {
        # in case of error, let user know what happend
        stop(e)
      }, finally = {
        # drop view
        idaDropView(tmpView)
      }
  )
  
  object.pred <- ida.data.frame(outtable)
  return(object.pred)
}

# ----------------------------------------------------------------------
