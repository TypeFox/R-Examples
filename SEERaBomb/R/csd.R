csd=function(seerSet,brkst=c(0),brksy=c(1973),brksa=c(0),trts=NULL,PYLong=FALSE,firstS="all"){ 
  #   surv=yrdx=modx=db=casenum=radiatn=cancer=trt=yrbrth=agedx=L2D=NULL
  # brkst=c(0);brksy=c(1973,2000);brksa=c(0);seerSet=pf
  print(binSt<-levels(cut(brkst+0.1,breaks=c(brkst,100)))) #this is just to make a vector of tsd interval/row names 
  print(binSy<-levels(cut(brksy,breaks=c(brksy,seerSet$yearEnd+1),right=FALSE,dig.lab=4))) # year at diagnosis break points
  print(binSa<-levels(cut(brksa+0.1,breaks=c(brksa,126)))) # age at diagnosis breaks
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
    for (R in trts) { 
      cat("treatment:",R,"\n")
      for (biny in binSy) {
        binIndxy=getBinInfo(biny,binSy)["index"]
        for (bina in binSa) {
          binIndxa=getBinInfo(bina,binSa)["index"]
          mids=NULL
          Obs=vector(mode="list",length=0)
          Exp=vector(mode="list",length=0)
          AgeO=vector(mode="list",length=0)  # mean age of observed cases
          PYA=vector(mode="list",length=0)
          if (PYLong) PYL=vector(mode="list",length=0)
          for (bint in binSt) {
            binIndxt=getBinInfo(bint,binSt)["index"]
            L1=post1PYOc(canc,brkst,binIndxt,brksy,binIndxy,brksa,binIndxa,Trt=R,PYLong=PYLong,yearEnd,firstS,secondS=secondS)
            Exp[[bint]]=getE(L1$LPYM,D,ageStart,ageEnd,yearEnd,firstS,secondS)
            Obs[[bint]]=L1$O[firstS,secondS,drop=FALSE]
            AgeO[[bint]]=L1$AgeO
            PYA[[bint]]=L1$PYA
            if (PYLong) PYL[[bint]]=L1$PYL
            mids=c(mids,L1$binMidPnt)
          } # loop on tsx bins
          SL[[R]][[biny]][[bina]]$mids=mids
          SL[[R]][[biny]][[bina]]$Obs=Obs
          SL[[R]][[biny]][[bina]]$Exp=Exp
          SL[[R]][[biny]][[bina]]$AgeO=AgeO
          SL[[R]][[biny]][[bina]]$PYA=PYA
          if (PYLong) SL[[R]][[biny]][[bina]]$PYL=PYL
        } # loop on agedx bins
      }# loop on yrdx bins
    } # loop on R
    SL
  })
  tsdn=paste0("b",paste(brkst,collapse="_"))
  seerSet$L[[tsdn]]=SL
  seerSet$active=N=length(seerSet$L) 
  seerSet$series=data.frame(index=1:N,series=names(seerSet$L))
  cat(paste("Current active series:",seerSet$active,"\n")) 
  print(seerSet$series) 
  seerSet$DF=getDF(seerSet)
  print(proc.time() - ptm)
  seerSet
}

