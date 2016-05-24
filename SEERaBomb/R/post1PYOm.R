post1PYOm=function(D,brks=c(0,2,5),binIndx=1,yearEnd) { survC=surv=agedx=py=year=ageL=NULL 
  # yearEnd=ceiling(max(D$yrdx+D$surv))
  # yearStart=floor(min(D$yrdx)) #no, stuck with 1973 as earliest year, else redesign fillPYM
  binS=levels(cut(brks+0.1,breaks=c(brks,100)))
  bin=binS[binIndx]
  LL=getBinInfo(bin,binS)["LL"]
  D=D%>%mutate(survC=cut(surv,breaks=c(-1,brks,100),include.lowest = TRUE)) 
  DO=D%>%filter(survC==bin) 
  O=sum(DO$status)
  D=D%>%filter(surv>LL) 
  D=D%>%mutate(py=getPY(surv,bin,binS,brks)) # getpy leaves zeros when surv end is left of LL
  D=D%>%filter(py>0)  #get rid of such rows upfront
  D=D%>%mutate(ageL=agedx+LL) 
  D$year=floor(D$yrdx+LL)
  D=D%>%select(py,ageL,year)
  binMidPnt=LL+sum(D$py)/dim(D)[1]/2
  PYin=as.matrix(D)
  yrs=1973:yearEnd; ages=0.5:125.5  # this + next line = initiate PYM with zeros
  PYM=matrix(0,ncol=length(yrs),nrow=length(ages),dimnames=list(ages,yrs)) 
  fillPYM(PYin,PYM)
  list(PYM=PYM,O=O,binMidPnt=binMidPnt)
}



