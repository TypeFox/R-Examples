makeData <- function(data = NULL, model.type = 1, modality = NULL, covariates = NULL){
  if(!is.data.frame(data)){
    stop("Data must be a data.frame!")
  }
  if(is.null(colnames(data))){
    stop("Column names of data are not given! 
         Please give the names to indicate \"TP\", \"FN\", \"FP\", \"TN\"!")
  }
  datanames = tolower(colnames(data))
  colnames(data) = datanames
  I = dim(data)[1]
  fic = c("tp","tn","fp","fn") # four important components
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  
  ###### check main data
  if(all(fic %in% datanames)){
    tf.integer = apply(cbind(data$tp,data$tn,data$fp,data$fn),2,function(x) all(is.wholenumber(x)))
    tf.pos = apply(cbind(data$tp,data$tn,data$fp,data$fn)>=0,2,function(x) all(x))
    if(all(tf.integer)){
      if(all(tf.pos)){
        if("studynames" %in% datanames){
          #message("Data is ok! Study Names are given!")
        }else{
          #print("Data is ok! Study names are not given!")
          data$studynames = paste("study[",c(1:I),"]",sep="")
        }
      }else{
        fvp = paste("data$",fic[tf.pos==FALSE],sep="",collapse=" and ")
        stop(fvp," has some negative value!!!")
      }
    }else{
      fvi = paste("data$",fic[tf.integer==FALSE],sep="",collapse=" and ")
      stop(fvi," has some non-integer value!!!")
    }
  }else{
    ffic = toupper(paste(fic[!(fic %in% datanames)],collapse=" "))
    stop(paste("Data is not ok!",ffic,"are missing!!!",sep=" "))
  }
  ########## check covariates
  if(is.null(covariates) || covariates==FALSE){
    cov.flag = FALSE
    #message("Covariates are not given!!!")
  }else{
    if(is.character(covariates) || is.numeric(covariates)){
      cov.flag = TRUE
      if(length(covariates)>=length(datanames)){
        cov.flag = FALSE
        stop("Covariates are given! Number of covariates can NOT be larger than number of variables!!! ")
      }
      if(is.numeric(covariates)){
        if(all(covariates %in% c(1:length(datanames)))){
          covariates = datanames[covariates]
        }else{
          cov.flag = FALSE
          stop("Covariates are given! Please give correct column number to indicate covariats!!!")
        }
      }
      if(!all(covariates %in% datanames)){
        cov.flag = FALSE
        stop("Covariates are given! Please give correct names of covariates!!!")
      }
      if(any(c(fic,"studynames") %in% covariates)){
        cov.flag = FALSE
        stop("Covariates are given! Covariates can not be \"studynames\", \"TP\", \"TN\", \"FP\", \"FN\"!!!")
      }
    }else{
      cov.flag = FALSE
      stop("Covariates are given! Argumens should either be character or integer!")
    } 
  }
  
  ########## check modality
  if(is.null(modality) || modality==FALSE){
    mod.flag = FALSE
  }else{
    if(is.character(modality) || is.numeric(modality)){
      mod.flag = TRUE
      if(length(modality)!=1){
        cov.flag = FALSE
        stop("Modality is given! Can Only have one modality !!! ")
      }
      if(is.numeric(modality)){
        if(modality %in% c(1:length(datanames))){
          modality = datanames[modality]
        }else{
          mod.flag = FALSE
          stop("Modality are given! Please give correct column number to indicate modality!!!")
        }
      }
      if(!all(modality %in% datanames)){
        mod.flag = FALSE
        stop("Modality are given! Please give correct names of modality!!!")
      }
      if(any(c(fic,"studynames") %in% modality)){
        mod.flag = FALSE
        stop("Modality are given! Mdality can not be \"studynames\", \"TP\", \"TN\", \"FP\", \"FN\"!!!")
      }
    }else{
      mod.flag = FALSE
      stop("Modality are given! Argumens should either be character or integer!")
    } 
  }
  if(any(modality %in% covariates)){
    stop("Modality can not be the same as covariates!!!")
  }
  
  if(mod.flag){
    um = as.character(unique(data[,modality]))
    ind = lapply(um, function(x){
      which(data[,modality]==x)
    })
    data_temp = lapply(ind, function(x){data[x,]})
    data = do.call(rbind, data_temp)
  }
  
  TP = data$tp
  TN = data$tn
  FP = data$fp
  FN = data$fn
  
  n1 = TP+FN
  n0 = FP+TN
  N = n1+n0
  
  Ntrials = matrix(rbind(n1,n0),2*I,1)
  studynames = rep(data$studynames,each=2)
  
  Y = matrix(rbind(TP,TN),2*I,1)*(model.type==1) +  ###### model.type 1 (se & sp)
    matrix(rbind(TP,FP),2*I,1)*(model.type==2) +    ###### model.type 2 (se & (1-sp))
    matrix(rbind(FN,TN),2*I,1)*(model.type==3) +    ###### model.type 3 ((1-se) & sp)
    matrix(rbind(FN,FP),2*I,1)*(model.type==4)      ###### model.type 4 ((1-se) & (1-sp))
  
  
  if(mod.flag){
    data.modality = data[,modality]
    um = as.character(unique(data.modality))
    munames = paste("mu.",um,sep="")
    nunames = paste("nu.",um,sep="")
    
    for(i in c(1:length(um))){
      ind = 1*(data.modality==um[i])
      mu.temp = matrix(rbind(ind,rep(0,I)),2*I,1)
      nu.temp = matrix(rbind(rep(0,I),ind),2*I,1)
      assign(munames[i], mu.temp)
      assign(nunames[i], nu.temp)
    }
    
    if(cov.flag){
      alphanames = paste("alpha.",covariates,sep="")
      betanames = paste("beta.",covariates,sep="")
      for(j in c(1:length(covariates))){
        cov.temp = data[[covariates[j]]]
        assign(covariates[j],cov.temp)
        alphatemp = matrix(rbind(cov.temp,rep(0,I)),2*I,1)
        betatemp = matrix(rbind(rep(0,I),cov.temp),2*I,1)
        assign(alphanames[j],alphatemp)
        assign(betanames[j],betatemp)
      } ##### end for loop
      modality.data = data.frame(studynames, Ntrials, Y, id=c(1:(2*I)), 
                                 mget(c(munames,nunames)),mget(c(alphanames,betanames)))
    }else{ ###### NO covariates
      modality.data = data.frame(studynames, Ntrials, Y, id=c(1:(2*I)), mget(c(munames,nunames)))
    }
    internaldata = modality.data
  } else{ ##### NO modality
    mu = matrix(rbind(rep(1,I),rep(0,I)),2*I,1)
    nu = matrix(rbind(rep(0,I),rep(1,I)),2*I,1)
    
    if(cov.flag){
      alphanames = paste("alpha.",covariates,sep="")
      betanames = paste("beta.",covariates,sep="")
      for(j in c(1:length(covariates))){
        cov.temp = data[[covariates[j]]]
        assign(covariates[j],cov.temp)
        alphatemp = matrix(rbind(cov.temp,rep(0,I)),2*I,1)
        betatemp = matrix(rbind(rep(0,I),cov.temp),2*I,1)
        assign(alphanames[j],alphatemp)
        assign(betanames[j],betatemp)
      } ##### end for loop
      covardata = data.frame(studynames, Ntrials, Y, id=c(1:(2*I)),
                             mu, nu, mget(c(alphanames,betanames)))
      internaldata = covardata
    }else{ ###### NO covariates
      internaldata = data.frame(studynames, Ntrials, Y, id=c(1:(2*I)), mu, nu)
    }
  }
  
  outdata = list()
  outdata$internaldata = internaldata
  outdata$originaldata = data
  outdata$covariates.setting = covariates 
  outdata$modality.setting = modality
  outdata$model.type = model.type
  
  return(outdata)   
}