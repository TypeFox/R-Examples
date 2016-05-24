tsd=function(seerSet,brks=c(0,2,5),trts=NULL,PYLong=FALSE,firstS="all"){ 
  #   surv=yrdx=modx=db=casenum=radiatn=cancer=trt=yrbrth=agedx=L2D=NULL
  print(binS<-levels(cut(brks+0.1,breaks=c(brks,100)))) #this is just to make a vector of tsd interval/row names 
  ptm <- proc.time()
  if(is.null(seerSet$L)) seerSet$L=list() # set list L to a blank list if it was never filled for this seerSet
  SL=with(seerSet, {  # subL (List within L), i.e. L is a list of time series sublists (SL). 
    #     yearEnd=max(popsa$year)
    SL=list() # initiate the SL that tsd() adds to L
    if (firstS[1]=="all") firstS=cancerS
    SL$firstS=firstS
    if (is.null(trts)) trts=levels(canc$trt)
    SL$trtS=trts
    print(trts)
    for (R in trts) 
    { 
      #     R="rad"
      print(R)
      mids=NULL
      Obs=vector(mode="list",length=0)
      Exp=vector(mode="list",length=0)
      AgeO=vector(mode="list",length=0)  # mean age of observed cases
#       AgeE=vector(mode="list",length=0)  # mean age of expected cases, i.e. of PY at risk
#       PY1=vector(mode="list",length=0)
      PYA=vector(mode="list",length=0)
#       PYT=vector(mode="list",length=0)
      if (PYLong) PYL=vector(mode="list",length=0)
      for (bin in binS) 
      {
        #       (bin=binS[1])
        #         print(bin)
        binIndx=getBinInfo(bin,binS)["index"]
        L1=post1PYO(canc,brks,binIndx,Trt=R,PYLong=PYLong,yearEnd,firstS,secondS=secondS)
        Exp[[bin]]=getE(L1$LPYM,D,ageStart,ageEnd,yearEnd,firstS,secondS)
        #          Obs[[bin]]=L1$O
        Obs[[bin]]=L1$O[firstS,secondS,drop=FALSE]
#         AgeE[[bin]]=L1$AgeE
        AgeO[[bin]]=L1$AgeO
#         PY1[[bin]]=L1$PY1
        PYA[[bin]]=L1$PYA
#         PYT[[bin]]=L1$PYT
        #         rws=rownames(L1$O)
        #         cols=colnames(L1$O)
        #         Obs[[bin]]=L1$O[rws%in%firstS,cols%in%secondS,drop=FALSE]
        if (PYLong) PYL[[bin]]=L1$PYL
        mids=c(mids,L1$binMidPnt)
      } # loop on tsx bins
      SL[[R]]$mids=mids
      SL[[R]]$Obs=Obs
      SL[[R]]$Exp=Exp
#       SL[[R]]$AgeE=AgeE
      SL[[R]]$AgeO=AgeO
      SL[[R]]$PYA=PYA
#       SL[[R]]$PY1=PY1
#       SL[[R]]$PYT=PYT
      if (PYLong) SL[[R]]$PYL=PYL
    } # loop on R
    SL
  })
  tsdn=paste0("b",paste(brks,collapse="_"))
  #   seerSet$bfn=paste0(seerSet$bfn,paste0("b",paste(brks,collapse="_")),txt)
  seerSet$L[[tsdn]]=SL
  seerSet$active=N=length(seerSet$L) 
  seerSet$series=data.frame(index=1:N,series=names(seerSet$L))
  cat(paste("Current active series:",seerSet$active,"\n")) 
  print(seerSet$series) 
  #   fn=paste0(outDir,"/",seerSet$bfn,".RData")
  #   cat("Writing L to ",fn,"\n")
  #   save(L,file=fn)
  #   seerSet$fL=fn
  print(proc.time() - ptm)
  seerSet
}

