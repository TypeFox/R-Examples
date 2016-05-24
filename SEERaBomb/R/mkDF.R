mkDF<-function(seerSet,srs=NULL) {
  age=trt=race=surv=year=py=popsa=cancer1=nO=O=E=ageE=ageEci=ageO=int=sig2O=ageOci=NULL 
  if (is.null(seerSet$L)) stop("seerSet L field is empty. Please run tsd on your seerSet object!") else {
    nms=names(seerSet$L)
    indx=seerSet$active
    cat("Current active series:",indx,"\n",sep="") 
    print(seerSet$series) 
    if (is.null(srs)) {indx=seerSet$active #grab most recent
                       cat(paste0("Using (active) series ",indx,", i.e. ", nms[indx],".\n")) 
                       #                        cat(paste0("Using active time series in L.   Index:",indx,"\n"))
    }
    if (length(srs)>1) {
      cat("... collapsing brks vector to a single string.\n") 
      srs=paste0("b",paste(srs,collapse="_"))
    }
    
    if (is.character(srs)){ indx=srs #grab ith most recent
                            cat(paste0("Using series ",which(indx==nms),", i.e. ", indx,".\n")) 
    }
    if (is.numeric(srs)){ indx=srs #grab ith most recent
                          cat("Using series ",indx,", i.e. ", nms[indx],".\n",sep="") 
    }
    
    L=seerSet$L[[indx]]
  }
  #   library(reshape2);library(plyr);library(SEERaBomb)
  #   load("~/Results/amlMDS/aveAgeFull.RData")
  #   seerSet=pm
  #   L=seerSet$L[[seerSet$active]]
  
  secondS=seerSet$secondS
  firstS=L$firstS
  #   options(stringsAsFactors=FALSE)
  f=data.frame(cancer1=firstS)
  #   f$cancer1
  fs=merge(f,data.frame(cancer2=secondS))%>%arrange(cancer1)
  #   fs$cancer1
  funf=function(x){
    y=left_join(f,x,by=c("cancer1"))
    y[is.na(y)]=0
    y
  }
  
  funfs=function(x,by){
    y=left_join(fs,x,by=by)
    y[is.na(y[,"n"]),"n"]=0
    y
  }
  
  funOE=function(x,by,out="O"){
    x=melt(x)
    names(x)=c("cancer1","cancer2",out)
    left_join(fs,x,by=by)
    #     y=left_join(fs,x,by=by)
    #     y[is.na(y[,"n"]),"n"]=0
    #     y
  }
  
  
  trts=L$trtS
  #   ae=NULL
  ao=NULL
  pya=NULL
  o=NULL
  e=NULL
  options(warn=-1)
  for (i in trts){
    #         i=trts[1]
    #     ae=rbind(ae,cbind(ldply(lapply(L[[i]]$AgeE,funfs,by=c("cancer1"))),trt=i)) 
    pya=rbind(pya,cbind(ldply(lapply(L[[i]]$PYA,funfs,by=c("cancer1"))),trt=i)) 
    ao=rbind(ao,cbind(ldply(lapply(L[[i]]$AgeO,funfs,by=c("cancer1","cancer2"))),trt=i)) 
    #     py=rbind(py,cbind(ldply(lapply(L[[i]]$PY1,funf)),trt=i)) 
    o=rbind(o,cbind(ldply(lapply(L[[i]]$O,funOE,by=c("cancer1","cancer2"),out="O")),trt=i)) 
    e=rbind(e,cbind(ldply(lapply(L[[i]]$E,funOE,by=c("cancer1","cancer2"),out="E")),trt=i)) 
  }
  # melt(L[[i]]$O[[1]])
  # lapply(L[[i]]$O,dim)
  # lapply(L[[i]]$E,dim)
  options(warn=0)
  #   head(ae,2)                       
  head(pya,2)                       
  head(ao,2)
  head(o,2)
  head(e,2)
  names(ao)[4:5]=c("ageO","sig2O")
  #   names(ae)[1]=c("int")
  names(pya)[1]=c("int")
  #   d=cbind(ae,ao%>%select(ageO:nO))
  #   d=ao%>%select(ageO:nO)
  #   py=as.data.frame(lapply(pya,rep,each=2))
  #   py=as.data.frame(lapply(py,rep,each=2))
  #   head(py,2)
  # d=cbind(d,py%>%select(py1:midPnt))
  #   d=cbind(d,py%>%select(py:t))
  d=cbind(pya%>%select(int:trt),ao%>%select(ageO:sig2O))
  d=cbind(d,O=o$O,E=e$E)
  d$int=factor(d$int,levels=unique(d$int))
  head(d)
  d=d%>%mutate(RR=O/E,
               rrL=qchisq(.025,2*O)/(2*E),
               rrU=qchisq(.975,2*O+2)/(2*E))
  #   glimpse(d)
  sem=sqrt(d$sig2/d$n)
  options(warn=-1) # expect warnings on t df of 0 and -1
  ci95=qt(0.975,d$n-1)*sem
  d$aeL=d$age-ci95
  d$aeU=d$age+ci95
  
  sem=sqrt(d$sig2O/d$O)
  ci95=qt(0.975,d$O-1)*sem
  options(warn=0) # turn off expect warnings on t df of 0 and -1
  d$aoL=d$ageO-ci95
  d$aoU=d$ageO+ci95
  i=sapply(d,is.numeric)
  d[i]=round(d[i],2)
  #   glimpse(d)
  d
} 
