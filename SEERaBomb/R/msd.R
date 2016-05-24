msd=function(canc,mrt,brks=c(0,2,5)){ #mortality since diagnosis (msd)
  surv=sex=O=E=NULL
  # yearEnd=ceiling(max(canc$yrdx+canc$surv))
  yearEnd=max(as.numeric(colnames(mrt$female)))
  canc=canc%>%filter(surv<200) # restrict to known survival times
  dm=canc%>%filter(sex=="male")%>%select(-sex)
  df=canc%>%filter(sex=="female")%>%select(-sex)
  d=list(male=dm,female=df)
  
  mrtF=mrt$female
  mrtF=rbind(mrtF,sapply(mrtF[111,],function(x) rep(x,15)))
  rownames(mrtF)=0:125
  mrtM=mrt$male
  mrtM=rbind(mrtM,sapply(mrtM[111,],function(x) rep(x,15)))
  rownames(mrtM)=0:125
  mrtM=mrtM[,as.character(1973:yearEnd)]
  mrtF=mrtF[,as.character(1973:yearEnd)]
  mrt=list(male=mrtM,female=mrtF)
  pts=c(male=dim(dm)[1],female=dim(df)[1],total=dim(dm)[1]+dim(df)[1])
  events=c(male=sum(dm$status),female=sum(df$status),total=sum(dm$status)+sum(df$status))
  Sexes=c("male","female")[pts[1:2]>0]
  print(Sexes)
  print(binS<-levels(cut(brks+0.1,breaks=c(brks,100)))) #this is just to make a vector of tsd interval/row names 
  DD=NULL
  for (S in Sexes) 
  { 
    print(S)
    mids=vector(mode="numeric",length=0)
    Obs=vector(mode="numeric",length=0)
    Exp=vector(mode="numeric",length=0)
    for (bin in binS) 
    {
      binIndx=getBinInfo(bin,binS)["index"]
      L1=post1PYOm(d[[S]],brks,binIndx,yearEnd)
      Exp[bin]=sum(L1$PYM*mrt[[S]]) 
      Obs[bin]=L1$O
      mids[bin]=L1$binMidPnt
    } # loop on tsx bins
    D=data.frame(int=factor(names(mids)),t=mids,O=Obs,E=Exp)
    D=D%>%mutate(RR=O/E,
                 rrL=qchisq(.025,2*O)/(2*E),
                 rrU=qchisq(.975,2*O+2)/(2*E),sex=S)
    DD=rbind(DD,D)
  } # loop on S (sexes)
  DD
}

