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

idaTree <- function(form, data, id, minsplit=50, maxdepth=10, qmeasure=NULL, 
    minimprove=0.01, eval=NULL, valtable = NULL, modelname=NULL) {
  
  modelName <- modelname
  
  # compute target variable varY from from
  
  ntab <- idaParseRFormula(form,data);
  varY  <- paste('\"',ntab$response,'\"',sep="")
  
  #Figure out whether this is regression or classification
  
  tableDef <- idaTableDef(data,F)
  
  targetType <- tableDef[tableDef$name==ntab$response,3];
  
  if(targetType=='NUMERIC') {
    isReg <- T;
  } else if(targetType=='CATEGORICAL') {
    isReg <- F;
  } else {
    stop('Type of target field not supported.')
  }
  
  # create view on required data dataTmp
  dataTmp <- data[,which(data@cols %in% c(ntab$cols, id,ntab$response))]
  
  model <- modelName
  
  if (is.null(model)) {
    model <- idaGetValidModelName('TREE_')
  } else {
    if(grepl(" ",model)) {
      stop("Space in model name not allowed.")
    }
    
    xx <- parseTableName(modelName);
    model <- paste('"',xx$schema,'"."',xx$table,'"',sep=''); 	
  }
  
  # check if given id is valid
  colu = dataTmp@cols
  if (!(id %in% colu))
    stop(simpleError(paste("Id variable is not available in ida.data.frame:", id)))
  
  id  <- paste('\"',id,'\"',sep="")
  
  valview <- NULL;
  tmpView <- NULL;
  tryCatch({	
        
        tmpView <- idaCreateView(dataTmp)
        
        if(!is.null(valtable)) {
          valview <- idaCreateView(valtable)
        }
        
        if(isReg) {
          
          callSP("IDAX.REGTREE ", model = model, intable = tmpView, id = id, target = varY, 
              minsplit = minsplit, maxdepth = maxdepth, qmeasure = qmeasure, 
              valtable = valview, minimprove = minimprove, eval = eval)
          
        } else {
          
          callSP("IDAX.DECTREE ", model = model, intable = tmpView, id = id, target = varY, 
              minsplit = minsplit, maxdepth = maxdepth, qmeasure = qmeasure, 
              valtable = valview, minimprove = minimprove, eval = eval)
        }					
        
      }, error = function(e, tmpView) {
        # in case of error, drop view and let user know, what happend
        stop(e)
      }, finally = {
        # drop views
        if(!is.null(tmpView))
          idaDropView(tmpView)
        
        if(!is.null(valview)) {
          valview <- idaDropView(valview)
        }
      }
  )
  
  idaRetrieveTreeModel(model);
}

idaRetrieveTreeModel <- function(modelName) {
  
  xx <- parseTableName(modelName);
  model <- xx$table
  modelSchema <- xx$schema
  
  modelTable <- paste('SELECT * FROM "',modelSchema,'"."',model,'_MODEL"',sep="")
  
  modelMain <- idaQuery(modelTable)
  
  isReg <- modelMain[1,1]=='regression'
  
  modelStr <- paste('SELECT * FROM "',modelSchema,'"."',model,'_NODES"',sep="")
  nodes <- idaQuery(modelStr)
  nodes$NODEID <- as.numeric(nodes$NODEID);
  nodes$PARENT <- as.numeric(nodes$PARENT);
  nodes$DEFAULTCHILD <- as.numeric(nodes$DEFAULTCHILD);
  
  modelStr <- paste('SELECT * FROM "',modelSchema,'"."',model,'_PREDICATES"',sep="")
  predicates <- idaQuery(modelStr)
  predicates$NODEID <- as.numeric(predicates$NODEID)
  
  modelStr <- paste('SELECT * FROM "',modelSchema,'"."',model,'_COLUMNS"',sep="")
  columns <- idaQuery(modelStr)
  
  if(!isReg) {
    modelStr <- paste('SELECT * FROM "',modelSchema,'"."',model,'_DISCRETE_STATISTICS"',sep="")
    discrete <- idaQuery(modelStr)
    discrete$NODEID <- as.numeric(discrete$NODEID)	
    discrete <- discrete[order(discrete[,"NODEID"]),]
  }
  
  # sort data by NODEID (no need to sort columns)
  nodes <- nodes[order(nodes[,"NODEID"]),]
  predicates <- predicates[order(predicates[,"NODEID"]),]
  
  # interprete CLASS and VALUE as factors, ISLEAF as logical
  
  if(!isReg) {
    nodes$CLASS <- as.factor(nodes$CLASS)
  } else {
    nodes$CLASS <- as.numeric(nodes$CLASS)
  }	
  
  if(!isReg) {
    discrete$VALUE <- as.factor(discrete$VALUE) # contains all possible classes
  }
  
  nodes$ISLEAF <- as.logical(nodes$ISLEAF)
  
  # compute availiable classes(=ylevels) and xlevels
  # create tree$terms
  
  if(!isReg) {
    classes <- levels(discrete$VALUE)
  } else {
    classes <- unique(nodes$CLASS)
  }
  
  target <- columns[columns$USAGETYPE=='predicted','COLUMNNAME'];
  contCols <- columns[columns$USAGETYPE=='active'&columns$OPTYPE=='continuous','COLUMNNAME']
  catCols <- columns[columns$USAGETYPE=='active'&columns$OPTYPE=='categorical','COLUMNNAME']
  
  allCols <- columns[(columns$USAGETYPE=='active')|(columns$USAGETYPE=='supplementary'),'COLUMNNAME']
  
  formStr <- paste(target,"~",paste(allCols,sep='',collapse='+'),sep='')
  form <- as.formula(formStr)
  
  pseudoData <- data.frame(xxxid=c(0,1));
  
  if(length(contCols)>0) {
    for(i in 1:length(contCols))
      pseudoData[contCols[i]]<-c(0,1)
  }
  
  if(length(catCols)>0) {
    for(i in 1:length(catCols))
      pseudoData[catCols[i]]<-c('a','b')
  }
  
  terms <- terms(form, data=pseudoData)
  
  xlevels <- attr(terms, "term.labels")
  
  # prepare tree$frame$var
  var <- matrix("", nrow=nrow(nodes))
  
  if(isReg) {
    levels(var) <- columns$COLUMNNAME
  }
  
  # compute tree$frame$n
  n <- as.integer(nodes$SIZE)
  
  # tree$frame$wt and tree$frame$dev are not given
  wt <- nodes$SIZE
  dev <- NA
  
  # Compute tree$frame$yval
  yval <- nodes$CLASS
  
  # tree$frame$complexity, tree$frame$ncompete and tree$frame$nsurrogate are not given
  complexity <- 0
  ncompete <- 0
  nsurrogate <- 0
  
  # compute tree$frame$yval2
  # sValue will be uses to compute tree$splits
  # compute dummy.xlevels in order to make tree print- and plottable
  
  if(!isReg) {
    yprob1 <- matrix(0, nrow=nrow(nodes), ncol=length(classes))
    yprob2 <- matrix(0, nrow=nrow(nodes), ncol=length(classes))
    yprob3 <- matrix(0, nrow=nrow(nodes))
    yprob4 <- matrix(0, nrow=nrow(nodes))
  }
  
  sValue <- matrix(0, nrow=nrow(nodes))
  # dummy.xlevels will contain relevant levels from every column
  dummy.xlevels <- list()
  # discretes contains the numbers of rows that are discrete
  discretes <- NULL
  
  for (i in 1:nrow(nodes)) {  	
    j <- predicates$NODEID[i]
    if(!nodes$ISLEAF[predicates$NODEID == j]) {
      var[i] <- predicates$COLUMNNAME[predicates$NODEID == 2*j]
      sValue[i] <- predicates$VALUE[predicates$NODEID == 2*j]
      if (predicates$OPERATOR[predicates$NODEID == 2*j] == "equal") {
        if (is.null(discretes)) {
          discretes <- i
        } else {
          discretes <- c(discretes, i)
        }
        if (!is.null(dummy.xlevels[[as.character(var[i])]])) {
          if (!(sValue[i] %in% dummy.xlevels[[as.character(var[i])]])) {
            dummy.xlevels[[as.character(var[i])]] <- c(dummy.xlevels[[as.character(var[i])]], sValue[i])
          }
        } else {
          dummy.xlevels[[as.character(var[i])]] <- c("<other>", sValue[i])
        }
      }
    } else {
      var[i] <- "<leaf>"
    }
    
    if(!isReg) {
      for (k in 1:length(classes)) {
        yprob1[i,k] <- discrete[discrete$NODEID==j & discrete$VALUE==classes[k],"COUNT"]
        yprob2[i,k] <- discrete[discrete$NODEID==j & discrete$VALUE==classes[k],"RELFREQUENCY"]
      }
      yprob3[i] <- as.numeric(yval[i])
      yprob4[i] <- nodes$RELSIZE[i]
    }
  }
  
  # if there are discrete variables, also create dummy.csplit in order to make tree print- and plottable
  if (!is.null(discretes)) {
    dummy.csplit <- matrix(2, nrow = length(discretes), ncol = max(sapply(dummy.xlevels, length)))
    for (l in 1:length(discretes)) {
      ll <- discretes[l]
      dummy.csplit[l,match(sValue[ll],dummy.xlevels[[var[ll]]])] <- 1
    }
    dummy.csplit[,1] <- 3
  } else {
    dummy.csplit <- NULL
  }
  
  var <- as.factor(var)
  
  if(!isReg) {		
    colnames(yprob4) <- "nodeprob"
    yval2 <- cbind(yprob3, yprob1, yprob2, yprob4)
  }	
  
  # create tree$frame
  
  if(isReg) {
    frame <- data.frame(var = var, n = n, wt = wt, dev = dev, yval = yval, complexity = complexity, ncompete = ncompete, nsurrogate = nsurrogate)
  } else {
    frame <- data.frame(var = var, n = n, wt = wt, dev = dev, yval = yval, complexity = complexity, ncompete = ncompete, nsurrogate = nsurrogate, yval2 = I(yval2))
  }
  
  
  rownames(frame) <- nodes$NODEID
  
  # put nodes into right order
  score <- (nodes$NODEID - 2^(floor(log(nodes$NODEID,2))  +0.5)) / 2^(floor(log(nodes$NODEID,2)+1))
  frame <- frame[order(score),]
  
  # compute tree$splits
  count <- nodes$SIZE
  ncat <- -1
  improve <- NA
  index <- suppressWarnings(as.numeric(sValue))
  adj <- NA
  
  splits <- cbind(count, ncat , improve, index, adj)
  rownames(splits) <- var
  
  ## if tree is a small tree (i.e. root node is leaf), don't cut the leaves
  smallTree = sum(nodes$ISLEAF == FALSE) == 1
  if (!smallTree) {
    splits <- splits[!nodes$ISLEAF,]
    splits <- splits[order(score[!nodes$ISLEAF]),]
    dummy.csplit <- dummy.csplit[order(score[discretes]),]
  }
  
  
  # add rpart functions print, summary, text to tree
  if(isReg) {
    dummy <- rpart(X~Y,data=data.frame(X=1,Y=1))$functions
  } else {
    dummy <-rpart(X~Y,data=data.frame(X="Q",Y=1))$functions
  }
  
  functions = list()
  functions$print <- dummy$print
  functions$summary <- dummy$summary
  functions$text <- dummy$text
  
  # create tree and set rpart-specific attributes
  tree <- list(frame = frame, terms = terms, splits = splits, functions = functions, model=modelName, modelTable=modelTable,method=ifelse(isReg,'anova','class'))
  attr(tree, "dummy.xlevels") <- dummy.xlevels
  attr(tree, "dummy.csplit") <- dummy.csplit
  
  # set attributes
  # format and dummies are required for proper printing / plotting
  attr(tree, "xlevels") <- xlevels
  attr(tree, "ylevels") <- classes
  
  attr(tree, "xlevels") <- attr(tree, "dummy.xlevels")
  tree$csplit <- attr(tree, "dummy.csplit")
  
  # adjust x$splits for discrete values
  # ncat <- 2 changes from continuous to discrete variable
  # index <- n selects row from x$csplit
  tree$splits[is.na(tree$splits[,"index"]),"ncat"] <- 2
  if (!is.null(tree$csplit)) {
    if (!is.matrix(tree$csplit)) {
      tree$csplit <- t(as.matrix(tree$csplit))
    }
    n.discrete <- nrow(tree$csplit)
    tree$splits[is.na(tree$splits[,"index"]),"index"] <- 1:n.discrete
  }
  
  # set class and return
  class(tree) <- c('idaTree','rpart')
  return(tree)
}

plot.idaTree <- function(x, ...) {
  
  if(x$method=='class') {
    prp(x, type = 2, extra = 104, nn = TRUE, varlen = 0, faclen = 0, shadow.col = "grey", fallen.leaves = TRUE, branch.lty = 3,...)
  } else {
    prp(x, type = 2, extra = 101, nn = TRUE, varlen = 0, faclen = 0, shadow.col = "grey", fallen.leaves = TRUE, branch.lty = 3,...)
  }
}

predict.idaTree <- function(object, newdata, id,...) {
  
  newData <- newdata
  
  outtable <- idaGetValidTableName(paste("PREDICT_",sep=""))
  
  colu = newData@cols
  if (!(id %in% colu))
    stop(simpleError(paste("Id variable is not available in ida.data.frame:", id)))
  
  id  <- paste('\"',id,'\"',sep="")
  tmpView <- idaCreateView(newData)
  
  tryCatch({	
        
        if(object$method=='class') {	
          callSP("IDAX.PREDICT_DECTREE ",
              model=object$model,
              intable=tmpView,
              id=id,
              outtable=outtable,
              ...)
        } else {
          callSP("IDAX.PREDICT_REGTREE ",
              model=object$model,
              intable=tmpView,
              id=id,
              outtable=outtable,
              ...)
        }
        
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
