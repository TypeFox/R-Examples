ndlCuesOutcomes <- function(formula, data, frequency=NA, numeric2discrete=function(x) Hmisc::cut2(x,g=g.numeric), g.numeric=2, check.values=TRUE, ignore.absent=NULL, variable.value.separator="", ...)
{

#  response = as.character(formula[2])
  response=gsub("[ ]+"," ",paste(deparse(formula[[2]],width.cutoff=500),collapse=""))
  predictors=gsub("[ ]+"," ",paste(deparse(formula[[3]],width.cutoff=500),collapse=""))
  n.predictors = length(predictors)

  Outcome = data[,response]

  if(predictors == ".")
    { predictors = colnames(data)
      predictors = predictors[predictors!=response]
      if(!is.na(frequency) & is.character(frequency) & length(frequency)==1)
        predictors = predictors[predictors!=frequency]
      formula = as.formula(paste(c(response,paste(predictors,collapse=" + ")),collapse=" ~ "))
    } 
  else
    { predictors = levels(as.factor(unlist(lapply(attr(terms.formula(formula),"term.labels"),function(x) strsplit(x,"[ ]*([\\+]|[\\|]|[:])[ ]*")))))
    }

  if(is.na(frequency[1]))
    Frequency = rep(1, length(Outcome))
  else
    { if(is.character(frequency) & length(frequency)==1)
        Frequency = data[,frequency]
      if(is.numeric(frequency) & length(frequency)==NROW(data))
        Frequency = frequency
    }

# Check for NA's

  na.predictors <- names(which(sapply(data[predictors], function(x) length(which(is.na(x))))!=0))
  if(length(na.predictors)!=0)
    stop(paste(c("The following predictors have NA's: ",paste(na.predictors, collapse=", ")), collapse=""))

# Discretize numeric predictors

  if(!is.null(numeric2discrete))
    { numeric.predictors <- names(which(sapply(data[predictors], is.numeric)))
      n.numeric.predictors <- length(numeric.predictors)
      if(n.numeric.predictors>=1)
        for(i in 1:n.numeric.predictors)
           data[numeric.predictors[i]] <- lapply(lapply(data[numeric.predictors[i]],numeric2discrete),function(x) gsub(" ","",x))
    }

# Prefix predictor names to predictor values

#  data <- apply(data[c(response,predictors)],c(1,2),as.character);
#  for(i in 1:n.predictors)
#     for(j in 1:nrow(data))
#        data[j,predictors[i]] <- paste(c(predictors[i],data[j,predictors[i]]),collapse=variable.value.separator)
#  data <- data.frame(data)

  data <- data.frame(cbind(data[response],
                     sapply(predictors,
                           function(i)
                                   sapply(data[,i],
                                         function(x) 
                                                 { if(check.values)
                                                     x = gsub("_",".",as.character(x))
                                                   else
                                                     if(length(grep("_",as.character(x)))!=0)
                                                       stop("Underscores '_' in predictor values: change manually or set 'check.values=TRUE'.")
                                                   if(is.null(ignore.absent))
                                                     paste(c(i,as.character(x)), collapse=variable.value.separator)
                                                   else
#                                                     ifelse(as.character(x)==ignore.absent,ignore.absent,paste(c(i,as.character(x)), collapse=variable.value.separator))
                                                     ifelse(as.character(x) %in% ignore.absent, ignore.absent[1], paste(c(i,as.character(x)), collapse=variable.value.separator))
                                                 }))))

# sanity check on unique level names for the predictors
#
#  valuelist = vector(mode="list")
#  cnt = 0
#  for (i in 1:ncol(data)) {
#    if (colnames(data)[i] %in% predictors) {
#      cnt = cnt + 1
#      vals = levels(factor(data[,i]))
#      vals = vals[vals != "NA"]
#      valuelist[[cnt]]=vals
#      names(valuelist)[[cnt]]=colnames(data)[i]
#    }
#  }
#  allValues = unlist(valuelist)
#  allValuesTab = table(allValues)
#  if (sum(allValuesTab==1) != length(allValues)) {
#    cat("There are duplicate levels for different predictors:\n")
#    for (i in 1:length(allValuesTab)) {
#      if (allValuesTab[i] > 1) {
#        cat(names(allValuesTab[i]), " ")
#      }
#    }
#    cat("\n")
#    cat("Adding column names to factor levels\n")
#
#    cnt = 0
#    for (i in 1:ncol(data)) {
#      if (colnames(data)[i] %in% predictors) {
#        cnt = cnt + 1
#        vals = paste(data[,i], colnames(data)[i], sep=".")
#        vals[grep("NIL", vals)]="NIL"
#        data[,i]=vals
#      }
#    }
#  }


  # the input for naive discriminative learning

#  colnms=colnames(data)
#  cnt = 0
#  for (i in 1:ncol(data)) {
#    if (colnms[i] %in% predictors) {
#      cnt = cnt + 1
#      if (cnt == 1) {
#        v = as.character(data[,colnms[i]])
#      } else {
#        v = paste(v, as.character(data[,colnms[i]]), sep="_")
#      }
#    }
#  }

  predictor.combinations <- lapply(attr(terms.formula(formula),"term.labels"), function(x) strsplit(x,":")[[1]])
  n.predictor.combinations <- length(predictor.combinations)

  v = NULL

  for(i in 1:nrow(data))
     { u = NULL;
       for(j in 1:n.predictor.combinations)
          if(length(predictor.combinations[[j]])>1)
             u = paste(c(u,paste(apply(data[i,predictor.combinations[[j]]],c(1,2),as.character),collapse=":")),collapse="_")
          else
             u = paste(c(u,as.character(data[i,predictor.combinations[[j]]])),collapse="_")
       if(!is.null(ignore.absent))
         { u = gsub(paste(c(":",ignore.absent[1],":"),collapse=""),":",u)
           u = gsub(paste(c("_",ignore.absent[1],":"),collapse=""),"_",u)
           u = gsub(paste(c(":",ignore.absent[1],"_"),collapse=""),"_",u)
           u = gsub(paste(c("^",ignore.absent[1],":"),collapse=""),"",u)
           u = gsub(paste(c(":",ignore.absent[1],"$"),collapse=""),"",u)
           u = paste(grep(ignore.absent[1],unique(strsplit(u,"_")[[1]]),value=TRUE,invert=TRUE),collapse="_")
         }
       v = c(v,u)
     }

#  for(i in 1:nrow(data))
#     v = c(v,paste(sapply(predictor.combinations,
#                         function(cols)
#                                 if(length(cols)>1)
#                                   paste(apply(data[i,cols],c(1,2),as.character),collapse=":")
#                                 else
#                                   as.character(data[i,cols])),
#           collapse="_"))

  if(!is.null(ignore.absent))  
    { n.all.cues.absent <- length(which(sapply(v, function(x) x=="")))
      if(n.all.cues.absent>0)
        stop(paste(c("All cues absent for N=", n.all.cues.absent, " lines in data."), collapse=""))
    }

  Cues = as.vector(v)

  cuesOutcomes = data.frame(Frequency=Frequency,
      Cues = Cues, Outcomes = as.character(data[,response]), 
      stringsAsFactors=FALSE)

#  if (maxn > 1) {
#    w=strsplit(cuesOutcomes$Cues, "_")
#    v = sapply(w, cueCoding.fnc, maxn)
#    cuesOutcomes$Cues = v
#  }

#  class(cuesOutcomes) <- "ndlCuesOutcomes"

  return(cuesOutcomes)

}
