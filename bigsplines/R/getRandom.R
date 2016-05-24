getRandom <- 
  function(formula,ndpts,data=NULL){
    
    # get random effect terms
    rpar <- parse(text=formula)
    if(length(rpar)>1) stop("Input 'random' must be one-sided formula.")
    rchar <- as.character(formula)[2]
    rterms <- unlist(strsplit(unlist(strsplit(rchar, "\\(?\\)", perl=TRUE)),"(?=\\()", perl=TRUE))
    ix <- which(rterms=="(")
    if(length(ix)>0) rterms <- rterms[-ix]
    ix <- which(rterms==" + ")
    if(length(ix)>0) rterms <- rterms[-ix]
    if(any(c(rterms==" - ",rterms==" * ",rterms==":",rterms=="/",rterms=="^"))) stop("Only additive random effects are allowed.")
    lenre <- length(rterms)
    
    # get random effect variables
    raneff <- reNames <- vector("list",lenre)
    for(j in 1:lenre){
      
      # split random effect terms
      spterms <- unlist(strsplit(rterms[j],"\\|"))
      if(length(spterms)!=2L) stop("Incorrect usage of grouping syntax in 'random' formula. See help files.")
      reNames[[j]] <- gsub(" ","",unlist(strsplit(spterms[1],"\\+")),fixed=TRUE)
      numre <- length(reNames[[j]])
      
      # get grouping variable
      if(j==1L){
        gname <- gsub(" ","",spterms[2],fixed=TRUE)
        grpvar <- data.frame(as.factor(eval(parse(text=gname),data,parent.frame())))
        if(nrow(grpvar)!=ndpts) stop("Length of variables in 'random' must match length of response.")
        names(grpvar) <- gname
      } else {
        gname <- gsub(" ","",spterms[2],fixed=TRUE)
        grpnew <- data.frame(as.factor(eval(parse(text=gname),data,parent.frame())))
        if(nrow(grpnew)!=ndpts) stop("Length of variables in 'random' must match length of response.")
        names(grpnew) <- gname
        grpvar <- cbind(grpvar,grpnew)
      }
      
      # get random effects
      if(numre==1L){
        
        if(reNames[[j]]!="1") {
          rename <- gsub(" ","",reNames[[j]],fixed=TRUE)
          raneff[[j]] <- data.frame(eval(parse(text=rename),data,parent.frame()))
          if(nrow(raneff[[j]])!=ndpts) stop("Length of variables in 'random' must match length of response.")
          names(raneff[[j]]) <- rename
        }
        
      } else {
        
        for(k in 1:numre){
          if(reNames[[j]][k]!="1") {
            rename <- gsub(" ","",reNames[[j]][k],fixed=TRUE)
            newraneff <- data.frame(eval(parse(text=rename),data,parent.frame()))
            if(is.factor(newraneff[,1])) stop("Incorrect usage of nesting syntax in 'random' formula. See help files.")
            if(nrow(newraneff)!=ndpts) stop("Length of variables in 'random' must match length of response.")
            names(newraneff) <- rename
            if(is.null(raneff[[j]])) { 
              raneff[[j]] <- newraneff 
            } else { 
              raneff[[j]] <- cbind(raneff[[j]], newraneff)
            }
          }
        } # end for(k in 1:numre)
        
      } # end if(numre==1L)
      
    } # end for(j in 1:lenre)
    
    return(list(grpvar=grpvar,reNames=reNames,raneff=raneff))
    
  }