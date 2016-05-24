  # This function proceeds Chi² tests crossing two sets of variables (Y and X) of the dataset
  "chisq.desc"=function(dataset,Y,X,method="proba",print=TRUE,report=FALSE,language="english"){
    assign("Y",Y,envir=.GlobalEnv)
    assign("X",X,envir=.GlobalEnv)
    c=names(dataset)[Y]
    d=names(dataset)[X]
    res=chisq.function(dataset,c,d,method,print)
    if(length(Y)<length(X)){
        a=X
        b=Y
    }else{
        a=Y
        b=X
    }
    a=a[!a%in%b]
    #f=names(dataset)[a]
    #g=names(dataset)[b]
    #e=c(a,b)
    tab_chi2(dataset,X,Y,0.05,print)
    if (report==TRUE){
        chisq.sweave(dataset,X,Y,language=language)
    }
  return(res)
  }


