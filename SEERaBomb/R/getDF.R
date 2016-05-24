getDF<-function(seerSet,srs=NULL) {
  age=trt=race=surv=year=py=popsa=cancer1=nO=O=E=ageE=ageEci=ageO=int=sig2O=ageOci=ageG=NULL 
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

  secondS=seerSet$secondS
  firstS=L$firstS
  #   options(stringsAsFactors=FALSE)
  f=data.frame(cancer1=firstS)
  fs=merge(f,data.frame(cancer2=secondS))%>%arrange(cancer1)

  funfs=function(x,by){
    y=left_join(fs,x,by=by)
    y[is.na(y[,"n"]),"n"]=0
    y
  }
  
  funOE=function(x,by,out="O"){
    x=melt(x)
    names(x)=c("cancer1","cancer2",out)
    left_join(fs,x,by=by)
  }
  
  trts=L$trtS
  ao=NULL
  pya=NULL
  o=NULL
  e=NULL
  options(warn=-1)
  for (i in trts){
    yrGrp=names(L[[i]])
    for (j in yrGrp) {
      ageGrp=names(L[[i]][[j]])
      for (k in ageGrp) {
        pya=rbind(pya,cbind(ldply(lapply(L[[i]][[j]][[k]]$PYA,funfs,by=c("cancer1"))),trt=i,yearG=j,ageG=k)) 
        ao=rbind(ao,cbind(ldply(lapply(L[[i]][[j]][[k]]$AgeO,funfs,by=c("cancer1","cancer2"))),trt=i)) 
        o=rbind(o,cbind(ldply(lapply(L[[i]][[j]][[k]]$O,funOE,by=c("cancer1","cancer2"),out="O")),trt=i)) 
        e=rbind(e,cbind(ldply(lapply(L[[i]][[j]][[k]]$E,funOE,by=c("cancer1","cancer2"),out="E")),trt=i)) 
      }
    }
  }
  options(warn=0)
  names(ao)[4:5]=c("ageO","sig2O")
  names(pya)[1]=c("int")
  d=cbind(pya%>%select(int:ageG),ao%>%select(ageO:sig2O))
  d=cbind(d,O=o$O,E=e$E)
  d$int=factor(d$int,levels=unique(d$int))
  d=d%>%mutate(RR=O/E,
               rrL=qchisq(.025,2*O)/(2*E),
               rrU=qchisq(.975,2*O+2)/(2*E))
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
  d=tbl_df(d)
  # print(glimpse(d))
  d
} 

