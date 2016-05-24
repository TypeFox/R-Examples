seerSet<-function(canc,popsa,Sex, Race="pooled",ageStart=15,ageEnd=85) {
  # gimic to get rid of unwanted notes in R CMD check
   agedx=age=age86=yrdx=sex=race=surv=modx=yrbrth=py=year=NULL 
  
  if (!"age"%in%names(popsa)) { #assume first time here
    # so rename age86 as age to help plotting and joining later
    popsa=popsa%>%mutate(age=age86)%>%select(-age86) 
  }
  
  if (!"age"%in%names(canc)) { #assume  first time here
    # let year be yrdx as a whole integer to free it to become a real
    canc=canc%>%mutate(year=yrdx) 
#when the next 2 mutations were in one call, yrdx would occasionally overwrite surv ... very weird.
    canc=canc%>%mutate(yrdx=round(yrdx+(modx-0.5)/12,3))    #modx=1=January 
    canc=canc%>%mutate(surv=round((surv+0.5)/12,3))%>%  
      select(-modx)%>%
      mutate(age=agedx+0.5) #convert ages at diagnosis to best guesses
    canc=canc%>%select(-agedx) 
  }  
  
#   canc=canc%>%filter(age>=(ageStart+0.5),age<(ageEnd+0.5),sex==Sex)
  canc=canc%>%filter(age>=ageStart,age<ageEnd,sex==Sex)
  popsa=popsa%>%filter(age>=ageStart,age<ageEnd,sex==Sex)
  if (Race!="pooled") {
    canc=canc%>%filter(race==Race)
    popsa=popsa%>%filter(race==Race) 
  }
  canc$cancer=factor(canc$cancer) # get rid of any opposite sex cancer type levels
  cancerS=levels(canc$cancer)
  
  if ("age86"%in%names(canc)) canc=canc%>%select(-age86) 
  if ("sex"%in%names(canc)) canc=canc%>%select(-sex) 
  if ("race"%in%names(canc)) canc=canc%>%select(-race) 
  if ("yrbrth"%in%names(canc)) canc=canc%>%select(-yrbrth)  
  # note: modx was removed from canc above
  
  if ("sex"%in%names(popsa)) popsa=popsa%>%select(-sex) 
  if ("race"%in%names(popsa)) popsa=popsa%>%select(-race)
  
  
  popsa=popsa%>%group_by(age,year)%>%summarize(py=sum(py))
  
  # and package it all up
  seerSet=list(canc=canc,popsa=popsa,sex=Sex,race=Race,ageStart=ageStart,ageEnd=ageEnd,cancerS=cancerS,yearEnd=max(popsa$year))
  class(seerSet)="seerSet"
  seerSet
} # return a list that can be attached or with-ed in other functions
