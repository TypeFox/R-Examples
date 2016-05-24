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

idaLm <- function(form, idadf, modelname = NULL, dropModel = TRUE, limit = 25){
  
  storeCovMat   <- TRUE
  clearExisting <- TRUE
              
  #SANITY-CHECKS
  if(!is.logical(storeCovMat)){
    stop("Parameter storeCovMat has to be logical.")
  }
  
  if(!is.logical(dropModel)){
    stop("Parameter dropModel has to be logical.")
  }
  
  if(dropModel&&(!is.null(modelname))) {
    stop("If dropModel is TRUE, you cannot set a model name.")
  }
  
  if(!is.logical(clearExisting)){
    stop("Parameter clearExisting has to be logical.")
  }
  
  if(!idaCheckProcedure("LINEAR_REGRESSION","idaLm",FALSE)){
    stop("Function not available.")                     
  } 
  
  if(!is.ida.data.frame(idadf)){
    stop("This function can only be applied to ida.data.frame objects.")
  }
  
  idaCheckConnection()
  formOut <- idaParseRFormula(form,idadf)
  if(length(formOut$response) == 0) {
    stop("No target variable specified.")
  }
  
  if(length(formOut$cols) > 78){
    stop("This function is currently not supported for more than 78 columns.")
  }
  
  #internal function that calls the stored procedure.
  idaLmSP <- function(formOut, idadf, modelname = NULL, storeCovMat = TRUE, dropModel = FALSE,
                      clearExisting = TRUE, ...){
    #CREATE MODELNAME
    if (is.null(modelname)) {
      modelname <- idaGetValidModelName('LINEAR_REGRESSION_')
    } else {    
      if(grepl(" ",modelname)) {
        stop("Space in modelname not allowed.")
      }   
      xx <- parseTableName(modelname);
      modelname <- paste('"',xx$schema,'"."',xx$table,'"',sep='')   
    }
      
    #ClearExisting
    if(idaModelExists(modelname)) { 
      if(clearExisting){
        idaDropModel(modelname)
      } else  {
        stop("Model already exists.")
      } 
    }

    #Create View
    input <- idaCreateView(idadf) 
    
    #Create formula with all attributes for the SP
    #Catching column-names with ' or " inside
    formOut$cols <- gsub('"', '""', formOut$cols)
    attr <- paste0("\"", formOut$cols, "\"", collapse=";")
    attr <- gsub("'", "''", attr)
    
    formOut$response <- gsub('"', '""', formOut$response)
    target <- paste0("\"", formOut$response, "\"")
    target <- gsub("'", "''", target)
    
    intercept <- "FALSE"
    if(formOut$intercept){
      intercept <- "TRUE"
    }
    #CALL STORED PROCEDURE
    tryCatch({ 
      res <- callSP("IDAX.LINEAR_REGRESSION ",
                    model = modelname,
                    intable = input,
                    target = target,
                    coldefrole = "ignore",
                    incolumn = attr,
                    storecovariancematrix = storeCovMat,
                    intercept = intercept)
      
    }, error = function(e) {
      # in case of error, let user know what happend
      stop(e)
    }, finally = {
      # drop view
      idaDropView(input)
    }
    )
    
    LinMod <- idaRetrieveidaLmModel(modelname)
      
    if(dropModel){
      idaDropModel(modelname)
    }
    return(LinMod)
  }
  
  idaLmSQL<-function(form,idadf, intercept){
    #calculates linear model for all continous variables and no in database model
    idaDotProductSQL <- function(idadf,colNames,intercept) {
    
      cols <- colnames(idadf);
      
      if(length(cols)>42)
        stop("Only up to 40 predictors are currently supported.");
      
      queryList <- c();	
      if(intercept) {	
        for(i in 1:length(cols)) {
          for(j in i:length(cols)) {
            queryList <- c(queryList,paste("SUM(\"", cols[i] ,"\"*\"", cols[j],"\") ",sep=''));
          }
          queryList <- c(queryList,paste("SUM(\"", cols[i] ,"\")",sep=''));
        }
        queryList <- c(queryList,"COUNT(*)");	
        
        query <- paste("SELECT ", paste(queryList,sep=',',collapse=',')," FROM ",idadf.from(idadf),sep='');
        res <- idaQuery(query);
        res[[length(res)]] <- as.numeric(res[[length(res)]])
        
        if (nrow(res[complete.cases(res),]) == 0) {
          stop("There is no valid data in the input table for building the model. Note that the rows with missing values in your input data are ignored.");
        }
        
        mdat <- matrix(1:(length(cols)+1)*(length(cols)+1),nrow=length(cols)+1,ncol=length(cols)+1,dimnames = list(c(colNames, "Intercept"),c(colNames, "Intercept")),byrow=T)
        
        r <- 1;
        c <- 1;
        for(i in 1:ncol(res)) {
          mdat[r,c] <- as.numeric(res[[i]][1]);
          mdat[c,r] <- mdat[r,c];
          c <- c + 1;
          if(c>length(cols)+1) {
            r <- r+1;
            c <- r;
          }
        }
        mdat
      } else {
        for(i in 1:length(cols)) {
          for(j in i:length(cols)) {
            queryList <- c(queryList,paste("SUM(\"", cols[i] ,"\"*\"", cols[j],"\") ",sep=''));
          }
        }
        
        query <- paste("SELECT ", paste(queryList,sep=',',collapse=',')," FROM ",idadf.from(idadf),sep='');
        res <- idaQuery(query);
        
        if (nrow(res[complete.cases(res),]) == 0) {
          stop("There is no valid data in the input table for building the model. Note that the rows with missing values in your input data are ignored.");
        }
        
        mdat <- matrix(1:(length(cols))*(length(cols)),nrow=length(cols),ncol=length(cols),dimnames = list(colNames,colNames),byrow=T)
        
        r <- 1;
        c <- 1;
        for(i in 1:ncol(res)) {
          mdat[r,c] <- as.numeric(res[[i]][1]);
          mdat[c,r] <- mdat[r,c];
          c <- c + 1;
          if(c>length(cols)) {
            r <- r+1;
            c <- r;
          }
        }
        mdat
      }  
    }
  
    
    #Creates a View with all columns that are in the formula form and deletes all
    #rows that have NA values. This might lead to an empty table!
    input <- idaCreateInputDataView(form,idadf)
    inputDf <- ida.data.frame(input$view)
    
    mat <- NULL
    tryCatch({
     mat <- idaDotProductSQL(inputDf,input$colNames,intercept)
     if(intercept){
       numrow <- mat[nrow(mat), ncol(mat)]
     }else{
       numrow <- nrow(inputDf)
     }
     newidaLm <- idaLmStatistics(mat, numrow = numrow)
    },error=function(e){
      stop(e)
    },finally={
      idaDropView(input$view)
    })
    
    card <- vector(mode = "integer", length = length(formOut$cols)) + 1
    names(card) <- formOut$cols
   
    newidaLm <- c(newidaLm, list(CovMat = mat, card = card, model = NULL))
    class(newidaLm) <- "idaLm"
    return(newidaLm)
  }

  #Checking for categorical columns so we know which method to choose.
  tableDef <- idaTableDef(idadf, FALSE)#all columns with their datatype
  tableDef <- tableDef[match(formOut$cols, tableDef$name), ]#look for the columns out of form
  numOfCatCols <- nrow(tableDef[tableDef$valType == "CATEGORICAL", ])
  
  if (numOfCatCols == 0 && dropModel && nrow(tableDef) < 42) {
   LinMod <- idaLmSQL(form, idadf, formOut$intercept)
  }else{
    LinMod <- idaLmSP(formOut, idadf, modelname = modelname, storeCovMat = storeCovMat,
                    dropModel = dropModel, clearExisting = clearExisting)
  }
  if(dropModel){
    LinMod$model <- NULL
  }
  return(LinMod)
}

idaLmStatistics<-function(mat, coeff = NULL, numrow){
  #This function calculates the statistics for the linear model.
  
  p <- ncol(mat)-1
  n <- numrow
  YTY <- mat[1, 1]
  XTX <- mat[-1, -1, drop = FALSE]
  XTY <- mat[1, -1]
  idaInvXTX <- NULL

  getLogLikelihood <- function(RSS, n) {
    -n/2*log(2*pi) - n/2*log(RSS/n) - n/2
  }
  
  getAIC <- function(L, n, p) {
    -2*L + 2*(p+1)
  }
  
  getBIC <- function(L, n, p) {
    -2*L + (p+1)*log(n)
  }
  
  getTest <- function(idaInvXTX, beta, n, p, sigmasq=1) {
    sdbeta = sqrt(diag(idaInvXTX)*sigmasq)
    tval   = beta/sdbeta
    pval   = 1 - abs(1 - 2*pt(tval, n-p))
    cbind(sdbeta, tval, pval)
  }
  
  getSigma <- function(RSS, df.residuals){
    sqrt(RSS / df.residuals)
  }
  
  #Try solve first, only if it fails use ginv
  try({idaInvXTX <- solve(XTX)}, silent = T)
  if(is.null(idaInvXTX)) {
    idaInvXTX <- ginv(XTX);
  }
  
  if(is.null(coeff)){
    coeff<-idaInvXTX %*% XTY
    rownames(coeff) = colnames(XTX)
  }
  beta <- coeff
  RSS  <- (YTY - 2* t(beta) %*% XTY + t(beta) %*% XTX %*% beta)[1, 1]
  names(RSS) <- NULL
  ran <- sum(abs(eigen(XTX, only.values = TRUE)$values) > 10^-8)#Every linear dependant attribute
  p    <- min(ran, n-1)                                         #has a zero eigenvalue.                   
  rank <- p
  df.residuals <- n-p
  likelihood <- getLogLikelihood(RSS, n)
  AIC   <- getAIC(likelihood, n, p)
  BIC   <- getBIC(likelihood, n, p)
  tests <- getTest(idaInvXTX, beta, n, p, RSS/(n-p))
  coefftab <- data.frame(Estimate = beta, Std.Error = tests[,1], t.value = tests[,2], p.value = tests[,3])      
  names(coefftab) <- c("Estimate", "Std.Error", "t.value", "Pr(>|t|)")
  rownames(coefftab)<-rownames(coeff)
  sigma <- getSigma(RSS, df.residuals)
  
  newidaLm <- list(coefficients = coeff, RSS = RSS, sigma = sigma, effects = NULL, rank = rank,
                   df.residuals = df.residuals, coefftab = coefftab, Loglike = likelihood,
                   AIC = AIC, BIC = BIC, numrow = numrow)
  return(newidaLm)
}

idaRetrieveidaLmModel<-function(modelname){
  #modelname with schema: "schema.table"
  #parse the modelname to get the table and schema names
  xx <- parseTableName(modelname)
  model <- xx$table
  modelSchema <- xx$schema
  
  
  #read all data from the model tables.
  modelquery <- paste0('SELECT * FROM "',modelSchema,'"."',model,'_MODEL"')
  Lmdata <- idaQuery(modelquery, na.strings = NULL)
  

  intercept <- 0
  interceptpos <- match(-1,Lmdata$VAR_ID)
  if(!is.na(interceptpos)){
    Lmdata$VAR_NAME[interceptpos] <- "Intercept" #changing it from (Intercept)
    intercept <- 1
  }
  coeff<-c()
  
  calcstats <- FALSE
  try({
    #Covariance-matrix
    covmatquery <- paste0('SELECT * FROM "',modelSchema,'"."',model,'_COVARIANCE_MATRIX"')
    covMat  <- idaQuery(covmatquery)
    nrowMat <- sqrt(nrow(covMat))#if full matrix saved
    #Change RCV-format into a matrix
    mat <- matrix(as.double(1:(nrowMat^2)),nrow=nrowMat,ncol=nrowMat)
    rowNames <- vector(mode = "character", length=nrowMat)
    covMat[, c(1, 3, 5)] <- data.matrix(covMat[, c(1, 3, 5)])
    for(i in 1:nrow(covMat)){
      mat[covMat$ROW_ID[i], covMat$COL_ID[i]] <- covMat$VAL[i]
      if(covMat$COL_ID[i] == 1){
        rowNames[covMat$ROW_ID[i]] <- covMat$ROW_NAME[i]
      }
    }
    colnames(mat) <- rowNames
    rownames(mat) <- rowNames
    
    #intercept is in line 2 instead of the last one.
    if(intercept){
      interceptpos <- match("Intercept", rowNames)
      switchIntercept <-  c((1:nrowMat)[-interceptpos], interceptpos)
      mat <- mat[switchIntercept, switchIntercept]
    }
    #orders the coeff table in the order of the CovMat if it exists
    coeff <- Lmdata$VALUE
    names(coeff) <- rownames(mat)[-1]#ASSUMPTION: Target variable in the first row
    
    if(!is.null(mat)){
      calcstats <- TRUE
    }
    
  },silent=TRUE)
  
  if(intercept){
    numrow <- mat[nrowMat, nrowMat]
  }else{
    modelinfo <- idaQuery(paste0("call idax.list_params('schema=", modelSchema, ", format=short,
                        where=modelname=''", model, "''')"))
    numrow <- as.numeric(modelinfo[match("valid_rows_in_intable", modelinfo$PARAMETERNAME), ]$PARAMETERVALUE)
  }
  
  #puts the coefficient of the model into a vector with the parsed names  
  for(i in 1:length(Lmdata$VAR_NAME)){
    if(!is.na(Lmdata$LEVEL_NAME[i])){
      tempname <- paste0(Lmdata$VAR_NAME[i],"___",Lmdata$LEVEL_NAME[i])
      coeff[tempname] <- Lmdata$VALUE[i]
    }else{
      coeff[Lmdata$VAR_NAME[i]] <- Lmdata$VALUE[i]
    }
  }
  coeff<-as.matrix(coeff)
  colnames(coeff) <- c("Coefficients")


  newidaLm <- list(coefficients = coeff, coefftab = NULL, RSS = NULL, sigma = NULL, effects = NULL,
                   rank = NULL, df.residuals = NULL, likelihood = NULL, AIC = NULL,
                   BIC = NULL, card = NULL, CovMat = NULL, model = modelname, numrow = numrow)
  
  if(calcstats){
    #Get the attribute names and their cardinality in order of the CovMat
    card <- table(Lmdata$VAR_NAME)[unique(Lmdata$VAR_NAME)]
    card <- card[!"Intercept" == names(card)]
   
    covmatnames <- rownames(mat)[-1]
    covmatnames <- covmatnames[!"Intercept" == covmatnames]
    
    #Get the original attribute names for all the dummy-variables
    splitted <- sapply(covmatnames, function(x){strsplit(x, split = "___")[[1]][1]})
    #find the first dummy-variable for each attribute
    startpos <- match(names(card), splitted)
    #reorder card in that order
    card <- card[order(startpos)]
    
    #Statistic functions from old idaLm function
    statistics <- idaLmStatistics(mat, coeff, numrow)
    #prepare object of the class idaLm
    newidaLm <- c(statistics, list(CovMat = mat, card = card, model = modelname))
  }
  class(newidaLm) <- c("idaLm")
  return(newidaLm)
}
  
print.idaLm <- function(x, ...) {
  cat("\nCoefficients:\n")
  if(is.null(x$coefftab)){
    print(x$coefficients)
  }else{
  coeff<-x$coefftab
  coeff$" " <-rep(" ", nrow(coeff))
  coeff[coeff$"Pr(>|t|)" < 0.1,   5] <- "."
  coeff[coeff$"Pr(>|t|)" < 0.05,  5] <- "*"
  coeff[coeff$"Pr(>|t|)" < 0.01,  5] <- "**"
  coeff[coeff$"Pr(>|t|)" < 0.001, 5] <- "***"
  print(coeff)
  cat("\n---")
  cat("\nSignif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
  cat("\nResidual standard error: ", x$sigma, " on ", x$df.residuals, " degrees of freedom\n\n")
  cat("Log-likelihood: ", x$Loglike, "\n")
  cat("AIC: ", x$AIC, "\n")
  cat("BIC: ", x$BIC, "\n")
  }
}

predict.idaLm <- function(object, newdata, id, outtable=NULL, ...){
  idaLm <- object
  #calculates the prediction of a given linear model for all datapoints in the input
  #
  #Args:
  #newdata:  ida.dat.frame that refers to a table containing data you want to be predicted
  #idaLm:    linear model  calculated by idaLm that has the same attributes as newdata
  #id:       name of the id column in newdata so we can match the prediction to the input
  #outtable: name of the table the results will be written in
  #
  #return:
  #ida.data.frame that refers to the table with the predicted values sorted by the id column
  if(is.null(outtable)){
    outtable <- idaGetValidTableName("PREDICT_LM_")
  }
  if(!"idaLm" %in% class(idaLm)){
    stop("object is not an object of class idaLm.")
  }
  if(!is.ida.data.frame(newdata)){
    stop("Parameter newdata is not an ida.data.frame.")
  }
  #newdata needs id column
  colu = newdata@cols
  if (!(id %in% colu))
    stop(simpleError(paste("Id variable is not available in ida.data.frame:", id)))
  
  #TODO: Add more verbose error message
  
  if(is.null(idaLm$model)){
    stop("Only models that are stored in-database can be applied. Re-train the model with the dropModel parameter set to FALSE to use the predict function.")
  }
  
  tmpview <- idaCreateView(newdata)
  #add quotation marks to support mixed case statements
  model    <- paste0("\"",idaLm$model,"\"")
  intable  <- paste0("\"",tmpview,"\"")
  id       <- paste0("\"",id,"\"")
  outtable <- paste0("\"",outtable,"\"")
  
  
  tryCatch({
    callSP("IDAX.PREDICT_LINEAR_REGRESSION",
           model    = model,
           id       = id,
           intable  = intable,
           outtable = outtable, ...)
    
  }, error = function(e) {
    # in case of error, let user know what happend
    stop(e)
  }, finally = {
    # drop view
    idaDropView(tmpview)
  }
  )
  object.pred <- ida.data.frame(outtable)
  return(object.pred)
}