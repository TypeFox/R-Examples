make.glm.formula <- function(formula, data, 
                                  name.runningtime = ".t",
                                  Min_T=Min_T,Max_T=Max_T, model=model, # used for formulasplit
                                ...){
# make formula for glm() on splitdata
# replace timevar by workingtime variable name

  special <- c("NPH","NLL", "NPHNLL") 
  Terms <- if (missing(data)){
    terms(formula, special, keep.order = TRUE)
  } else {
    terms(formula, special, data = data, keep.order = TRUE)
  }


  NamesNPHVars<- all_specials_vars( Terms,
                                   specials=c("NPH"),
                                   unique = TRUE,
                                   order="formula")
  
  NamesNPHNLLVars<- all_specials_vars( Terms,
                                      specials="NPHNLL",
                                      unique = TRUE,
                                      order="formula")
  

  

  modified <- 0
  newtermlabels <- labels(Terms)


      # change timevar in NPH() call
  if(length(NamesNPHVars) >0){
    for (i in attr(Terms, "specials")["NPH"]){
      for (k in 1:length(i)){        
        thecall <-  match.call(NPH, attr(Terms,"variables")[[i[k]+1]])
        modified <- modified + 1
        thecall[["timevar"]] <- as.name(name.runningtime)
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[k]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
      }
    }
  }
  
  
      # change timevar in NPHNLL() call
  if(length(NamesNPHNLLVars) >0){
    for (i in attr(Terms, "specials")["NPHNLL"]){
      for (k in 1:length(i)){        
        thecall <-  match.call(NPHNLL, attr(Terms,"variables")[[i[k]+1]])
        modified <- modified + 1
        thecall[["timevar"]] <- as.name(name.runningtime)
        indxterm <- variable2term(i, Terms)
        charcall<-deparse(thecall, 500)
        oldtermlabel <- newtermlabels[indxterm[k]]
        newtermlabels <- gsub(oldtermlabel, charcall, newtermlabels, fixed=TRUE)
      }
    }
  }
  
  
  if(modified > 0){
    formula <- reformulate(newtermlabels,
                           response = if (attr(Terms, "response")){ 
                             Terms[[2L]]
                           }
                           else NULL,
                           intercept = attr(Terms, "intercept"))
  }
  
  
  # Add missing default arguments
  #---------------------------------------------------------------------------------------- 
  fbis <- attr(terms(formula, keep.order = TRUE),"term.labels")

  if (model=="multiplicative" &  (length(attr(Terms, "specials")$NPHNLL)!=0)) {
    for (i in 1:(length(fbis))) {
      
      if ((substr(fbis[i],1,4)=="NPH(") | (substr(fbis[i],1,4)=="NLL(")) {
        fbis[i] <- paste((substring(fbis[i],first=1,last=(nchar(fbis[i]))-1)),
                         ")", sep="")
      }        
      
      if ((substr(fbis[i],1,7)=="NPHNLL(")) {                
        fbis[i] <- paste((substring(fbis[i],first=1,last=(nchar(fbis[i]))-1)),
                         ")",sep="")
        
        if (length(strsplit(fbis[i],"Boundary.knots")[[1]])==length(strsplit(fbis[i],"Boundary.knots.t")[[1]]) &
              length(strsplit(fbis[i],"Boundary.knots.t")[[1]])==1) {
          fbis[i] <- paste(substring(fbis[i],first=1,last=(nchar(fbis[i]))-1),
                           ",Boundary.knots.t=range(c(",Min_T,",",Max_T,")))",sep="")
        }                
        if (length(strsplit(fbis[i],"Boundary.knots.t")[[1]])==1) {
          fbis[i] <- paste(substring(fbis[i],first=1,last=(nchar(fbis[i]))-1),
                           ",Boundary.knots.t=range(c(",Min_T,",",Max_T,")))",sep="")
        }
        if (length(strsplit(fbis[i],"Boundary.knots")[[1]])==2) {
          fbis[i] <- paste(substring(fbis[i],first=1,last=(nchar(fbis[i]))-1),",Boundary.knots=range(",
                           strsplit(strsplit(as.character(fbis[i]),"(",fixed=TRUE)[[1]][2],",")[[1]][1],
                           "))",sep="")
        }  
      }
    }
  } 
  
return(formula)
  
}
