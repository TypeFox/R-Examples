mk2D<-function(seerSet, knots=5, write=FALSE, outDir="~/Results", txt=NULL,secondS=NULL) {
  if(!file.exists(outDir))  {   print(paste("Creating directory",outDir))
                                dir.create(outDir,recursive=TRUE)    }
  #   require(dplyr)
  #   require(mgcv) # Mixed GAM Computation Vehicle with GCV/AIC/REML smoothness estimation
  ptm <- proc.time()
  seerSet=with(seerSet, {
    if (is.null(secondS)) secondS=cancerS
    seerSet$secondS=secondS
    L2D=vector(mode="list",length=length(secondS)) 
    names(L2D)=secondS
    L2Dp=L2D; D=NULL; 
    knotsIn=knots
    for (i in secondS) {
      knots=knotsIn
      if (i=="ALL" & ageStart<15) knots=20   # make some knot numbers conditional, like this
      if (i=="otherCIS" & knots<10) knots=10   # 5 doesn't cut it for this, 10 does
      if (i=="liver" & knots<10) knots=10   # same here, to capture the cohort effect ripple
      if (i=="APL" & knots>5) knots=5   # here too many knots dive down too much into early calendary years of too few cases
      #the next if else deals with MDS and CMML starting only later. Not sure what other cancers are like this.
# if (i=="MDS")  {d=canc%>%filter(cancer%in%i,year>2000)
      if (i%in%c("MDS","MPN","RARS","AMLti","unknown"))  {d=canc%>%filter(cancer%in%i,year>2000)
                                 ps=popsa%>%filter(year>2000)
      } else 
        if (i=="CMML")  {d=canc%>%filter(cancer%in%i,year>1985) 
                         ps=popsa%>%filter(year>1985)
        } else  {
          d=canc%>%filter(cancer%in%i)
          ps=popsa
        }
      d=d%>%mutate(age=floor(age)+0.5)
      d=d%>%group_by(year,age)%>%summarise(cases=n())
      ps=ps%>%group_by(year,age)%>%summarise(py=sum(py)) 
      X=left_join(ps,d)%>%mutate(incid=1e5*cases/py)
      X[is.na(X)]=0
      cat(knots," knots, working on cancer ",i,":\n")
      L2D[[i]]=gam(cases ~ s(age)+s(year)+ti(age,year,k=knots)+offset(log(py)),
                   family=poisson(),data=X,method="REML") 
      prd=as.numeric(predict(L2D[[i]]))
      L2Dp[[i]]=cbind(cancer=i,X,Ecases=exp(prd),Eincid=1e5*exp(prd)/X$py  )
      D=rbind(D,L2Dp[[i]])
    } # i loop on cancers in secondS
    bfn<-paste0(substr(race,1,1),toupper(substr(sex,1,1)),"s",ageStart,"e",ageEnd,txt) #base of file names
    seerSet$D=tbl_df(D)
    seerSet$bfn=bfn
    if (write) {
      fL2D<-paste0(outDir,"/",bfn,"L2D.RData");    
      cat("writing L2D file ...",fL2D,":\n")
      save(L2D,file=fL2D)  
      seerSet$fL2D=fL2D
    }
    seerSet
  }) # end with
  print(proc.time() - ptm) 
  seerSet # return extended seerSet, now including D and file base name, and the L2D file name, if written
} #end func
