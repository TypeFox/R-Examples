Lexispoints <-
function (Bdata,transition,title,cov,legend)
{ z<- check.par (Bdata)
 # require (Epi)
  if (missing(cov)) cov=NA
  if (missing(title)) title <- "Title missing" 
  if (missing(legend)) legend <- "topleft"
 # year <- data.frame(YearTrans(Bdata))
  Bdata2 <- date_b(Bdata=Bdata,selectday=1,format.out="year")  
  zx <- TransitionAB(Bdata2,transition) # zx: nrow = number of transition selected
  Bdata2 <- subset(Bdata2,!is.na(zx$year))
  ID <- zx$id
  entry <- subset(Bdata2$start,Bdata2$ID%in%ID)
  exit <-  subset (zx$year,!is.na(zx$year)) # event of interest
  birth <- subset(Bdata2$born,Bdata2$ID%in%ID)
  istrans <- rep(1,length(ID))
  #zbb <- subset(Bdata2,!is.na(zx$year))
   Lcoh1 <- Lexis( id = ID,
               entry = list( CalTime=entry ),
               exit  = list( CalTime=exit, Age=exit-birth ),
               exit.status = istrans,
               data=Bdata2[Bdata2$ID%in%ID,], # zbb,
               merge=TRUE)       # 1900 added for x-axis
               # lex.dur = exit - entry
               # CalTime =  Lcoh1[,1] = date at entry (1900+)
               # Age = age at ENTRY
  # if (date_in_month_or_year == 1) moye <- 12 else moye <- 1
  moye <- 1
  AgeLow1 <- 5*trunc(min(na.omit(exit-birth)/(5*moye))) 
  AgeHigh1 <- 5*(trunc(max(na.omit(exit-birth)/(5*moye)))+1)
  PerLow1 <- 5* trunc(min(na.omit(exit)/(5*moye)))
  PerHigh1 <- 5* (trunc(max(na.omit(exit)/(5*moye)))+1)
  PerHigh1[PerHigh1-PerLow1 < AgeHigh1-AgeLow1] <- PerLow1 + AgeHigh1 - AgeLow1
  AgeHigh1[AgeHigh1-AgeLow1 < PerHigh1-PerLow1] <- AgeLow1 + PerHigh1 - PerLow1
  Lexis.diagram (age=c(AgeLow1,AgeHigh1),alab="Age",date=c(PerLow1,PerHigh1),dlab="Calendar Time",
  int=5,lab.int=5,col.life="black",lwd.grid=1,las=1,main=title)
  color <- palette (rainbow(7))
  Lcoh1$one <- 1
  if (is.na(cov)) 
   { poscov <- which (colnames(Lcoh1)%in%"one")
     covcat <- 1 } else
   { poscov <- which (colnames(Lcoh1)%in%cov) 
     covcat <- sort(unique(Lcoh1[,poscov]))
   }
  for (i in 1:length(covcat))
  { points(Lcoh1$CalTime[Lcoh1[,poscov]==covcat[i]]+Lcoh1$lex.dur[Lcoh1[,poscov]==covcat[i]],
     Lcoh1$Age[Lcoh1[,poscov]==covcat[i]]+Lcoh1$lex.dur[Lcoh1[,poscov]==covcat[i]], 
     pch=c(NA,16)[Lcoh1$lex.Xst[Lcoh1[,poscov]==covcat[i]]+1], col=color[i], cex=0.5 )
  }
 legend(legend, title=cov,legend=covcat[1:length(covcat)],pch=16,col=color[1:length(covcat)],bg="white")
 return (Lcoh=Lcoh1)
}
