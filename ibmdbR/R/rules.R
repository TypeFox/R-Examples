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

idaArule <- function(data,tid,item,maxlen=5,maxheadlen=1,minsupport=NULL,minconf=0.5,nametable=NULL,namecol=NULL, modelname=NULL) {
  
  if(!idaCheckProcedure("ASSOCRULES","idaArule",F)) {
    stop("idaArule is not available for the current connection.")
  }
  
  model <- modelname
  
  if (is.null(model)) {
    model <- idaGetValidModelName('ARULE_')
  } else {
    if(grepl(" ",model)) {
      stop("Space in model name not allowed.")
    }
    
    xx <- parseTableName(modelname);
    model <- paste('"',xx$schema,'"."',xx$table,'"',sep=''); 
    
    if(idaModelExists(model)){
      stop('Model with the same name already exists.')
    }
  }
  
    # create view from given data argument
    tmpView <- idaCreateView(data)
    
    tryCatch({		
          # validate given id and item argument
          colu = data@cols
          if (!(tid %in% colu)) {
            stop(simpleError(paste("tid variable is not available in ida.data.frame:", tid)))
          }
          if (!(item %in% colu)) {
            stop(simpleError(paste("item variable is not available in ida.data.frame:", item)))
          }
          
          tid <- paste('\"', tid, '\"', sep="")
          item <- paste('\"', item, '\"', sep="")
          
          # call stored procedure
          callSP("IDAX.ASSOCRULES",
              model=model,
              intable=tmpView,
              tid=tid,
              item=item,
              maxlen=maxlen,
              maxheadlen=maxheadlen,
              minsupport=minsupport,
              minconf=minconf,
              nametable=nametable,
              namecol=namecol)
        }, error = function(e) {
          # in case of error, let user know what happend
          stop(e)
        }, finally = {
          # drop view
          idaDropView(tmpView)
        }
    )
    
    result <- idaRetrieveRulesModel(model);
    
    return(result)
  }
   

idaRetrieveRulesModel <- function(modelname) {

  xx <- parseTableName(modelname);
  model <- xx$table
  modelSchema <- xx$schema
  
  itemTable <- paste('SELECT * FROM "',modelSchema,'"."',model,'_ITEMS" ORDER BY ITEMID',sep="")
  itemsetTable <- paste('SELECT * FROM "',modelSchema,'"."',model,'_ASSOCPATTERNS"',sep="")
  
  ruleTable <- paste('SELECT A.*,B.SUPPORT FROM "',modelSchema,'"."',model,'_ASSOCRULES" A',',"',modelSchema,'"."',model,'_ASSOCPATTERNS_STATISTICS" B WHERE A.ITEMSETID=B.ITEMSETID ORDER BY RULEID',sep="")
  headSupport <- paste('SELECT A.RULEID,B.SUPPORT FROM "',modelSchema,'"."',model,'_ASSOCRULES" A',',"',modelSchema,'"."',model,'_ASSOCPATTERNS_STATISTICS" B WHERE A.HEADID=B.ITEMSETID ORDER BY RULEID',sep="")
  
  item.out <- idaQuery(itemTable)
  itemset.out <- idaQuery(itemsetTable)
  rule.out <- idaQuery(ruleTable)
  head.out <- idaQuery(headSupport)
  
  rule.out$LIFT <- rule.out$CONFIDENCE/head.out$SUPPORT
  
  # how many itemsets are there?
  n.sets <- nrow(itemset.out)
  
  # itemsets are saved in sparse matrices, i. e. for each itemset
  # there is one row, containing one column for each item, if this item
  # is inside the itemset by true or false
  # need to convert itemids into position ids in order to create the sparse matrices
  # position[ ,1] contains the itemset id
  # position[ ,2] contains the item id
  position <- itemset.out[,c('ITEMSETID','ITEMID')]
  
# create sparse matrix with all itemsets
  sets <- sparseMatrix(position[,1],position[,2])
# load itemnames (only required ones 1:ncol(sets))
  sets.labels <- item.out[1:ncol(sets),3]
  
  if(length(unique(sets.labels))!=length(sets.labels)) {
    warning("Duplicates in item names are not supported. Using item ids.")
    sets.labels <- item.out[1:ncol(sets),2]
  }
  
# for arules compatibility, make sure labels are in 'AsIs' character format
  sets.labels <- I(as.character(sets.labels))
  
# create quality table
  quality <- rule.out[ ,c("SUPPORT", "CONFIDENCE", "LIFT")]
# for arules compatibility, columnnames must be lowercase
  names(quality) <- tolower(names(quality))
  
 info <- data.frame(data="RETREIVED FROM EXISTING MODEL", NA, support=0, confidence=0, model=I(model))
  
# get left and right hand sides from sets
# NOTE: calling as-Matrix again is necessary when there is only one rule
  if (nrow(rule.out)==1) {
    lhs.data <- as(Matrix(sets[rule.out[ ,"BODYID"], ]),"nsparseMatrix")
    rhs.data <- as(Matrix(sets[rule.out[ ,"HEADID"], ]),"nsparseMatrix")
  } else {
    lhs.data <- Matrix::t(sets[rule.out[ ,"BODYID"], ])
    rhs.data <- Matrix::t(sets[rule.out[ ,"HEADID"], ])
  }
# create lhs and rhs as "itemMatrix" objects, transpose lhs and rhs for compatibility
  lhs <- new("itemMatrix", data=lhs.data, itemInfo=data.frame(labels=sets.labels))
  rhs <- new("itemMatrix", data=rhs.data, itemInfo=data.frame(labels=sets.labels))
  
# create "rules" object
  rules <- new("rules", lhs=lhs, rhs=rhs, quality=quality, info=info)
  
  return(rules)
}

idaApplyRules <- function(modelname, newdata, tid, item) {
  
  outtable <- idaGetValidTableName(paste("APPLYRULES_",sep=""))
  
  colu = newdata@cols
  if (!(item%in% colu))
    stop(simpleError(paste("Item variable is not available in ida.data.frame:", item)))
  
  if (!(tid %in% colu))
    stop(simpleError(paste("TID variable is not available in ida.data.frame:", tid)))
  
  
  item  <- paste('\"',item,'\"',sep="")
  tid  <- paste('\"',tid,'\"',sep="")
  
  tmpView <- idaCreateView(newdata)
  
  tryCatch({	
        callSP("IDAX.PREDICT_ASSOCRULES",
            model=modelname,
            intable=tmpView,
            item=item,
            tid=tid,
            outtable=outtable)
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



