Lexislines.episodes <-
function (Bdata,Dlong,subjectsID,title)
{ # From lifelines.r
  # subjectsID <- Bdata$ID
  #  subjectsID <- c(1,19,20,208)
 # require (Epi)
  # Convert dates in years
  if (missing(subjectsID)) subjectsID <- sample(Bdata$ID,5,replace=FALSE)
  if (missing(title)) title <- "Title missing"
   if (is.null(attr(Bdata,"format.date"))) {print ("Lexisines.episodes: format.date is missing (attribute of data)")}
  z<- check.par (Bdata)
  if (missing(Dlong)) 
       { print ("Getting data in long file format. Patience please.",quote=FALSE)
       	 Dlong <- Biograph.long(Bdata)
       	 Dlong2 <- Dlong$Depisode
       	 print ("Data in long format produced. Lexis continues.") } else Dlong2 <- Dlong # on input, Dlong is Dlong#Depisode
       	 
  if (!is.data.frame(Dlong2)) {stop ("Dlong$Depisode is not a data frame. Please check") }
  namstates <- attr(Dlong2,"param")$namstates
   format.in <- attr(Dlong2,"format.date")
   format.born <- attr(Dlong2,"format.born")
   y <- date_convert (Dlong2$Tstart,format.in=format.in,format.out="year",format.born=format.born)
   Dlong2$TstartY <- y
   y <- date_convert (Dlong2$Tstop,format.in=format.in,format.out="year",format.born=format.born)
   Dlong2$TstopY <- y
   y <-date_convert (Dlong2$born,format.in=format.in,format.out="year",format.born=format.born) 
   bt <- y
   Dlong2$Tstartage <- Dlong2$TstartY - bt
   Dlong2$Tstopage <- Dlong2$TstopY - bt
         
  en1 <- Dlong2$TstartY
  ex1 <-  Dlong2$TstopY
  # Check whether duration is non-negative
  duration <- ex1 - en1
  duration.neg <- length (duration[duration < 0]) # number of negative durations
  if (duration.neg > 0)
  { print ("Lexislines.episodes.R: some durations are negative.")
    print (Dlong2[duration<0,])
    return
  }
  Lcoh1 <- Lexis( id = Dlong2$ID,
               entry = list( CalTime=en1 ),
               exit  = list( CalTime=ex1, Age=ex1-bt ),
               exit.status = Dlong2$status,
               data=Dlong2,
               merge=TRUE)
  nyear <- 5 
  if (max(na.omit(ex1-bt))-min(na.omit(en1-bt))<20) nyear <- 1  
  AgeLow <- nyear*trunc(min(na.omit(en1-bt))/nyear)
  AgeHigh <- nyear*trunc(max(na.omit(ex1-bt))/nyear+1)
  PerLow <- nyear* trunc(min(na.omit(en1)/nyear))
  PerHigh <- nyear* (trunc(max(na.omit(ex1)/nyear))+1)
  PerHigh[PerHigh-PerLow < AgeHigh-AgeLow] <- PerLow + AgeHigh - AgeLow
  AgeHigh[AgeHigh-AgeLow < PerHigh-PerLow] <- AgeLow + PerHigh - PerLow

  subjectsID2 <- subjectsID[subjectsID %in% as.numeric(Dlong2$ID)]  # delete IDs that do not exist
#  print (c(subjectsID2,Lcoh1$lex.id[1:100]))
#str(subjectsID)
#str(subjectsID2)
#str(Lcoh1$lex.id)
  
  #subjectsID <- subjectsID[Bdata$ns[Dlong2$ID %in%zz] %in% c(2:20)]  # delete IDs without transitions
  #Lcoh11 <- subset(Lcoh1, subjectsID2 %in% lex.id)
  Lcoh11 <- Lcoh1[Lcoh1$ID %in% subjectsID2,]
  #print (subjectsID2)
  # class (subjectsID2)
  
  colours <- rainbow(length(namstates))
  # plot lifelines
  Lcoh11$col<- 1
  for (i in 1:nrow(Lcoh11))
  { Lcoh11$col[i] <- colours[grep(Lcoh11$OR[i],namstates)]  }
  plot.Lexis( Lcoh11, grid=0:20*nyear, col=Lcoh11$col, xlim=c(PerLow,PerHigh),
     ylim=c(AgeLow,AgeHigh), lwd=2, las=1,col.grid="gray",
     main=title)
  
 # Mark the location of the events in the Lexis diagram
  pchh <- c(19,16,18)
  #Dlong22 <- subset(Dlong2,Dlong2$ID %in%subjectsID)
# ------------------------------------------------
# To prevent NOTE: no visible binding for global variable
# http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check
# See also Lexis.lines
TstopY <- Tstopage <- DES <- lex.id <- NULL
# -------------------------------------------------
  if (length(subjectsID2) < 10) points( Lcoh11$TstopY,Lcoh11$Tstopage,pch=substr(Lcoh11$DES,1,1), cex=0.7 )
 # Display ID 
  if (length(subjectsID2) < 20) 
   {  # Lcoh12 <- subset (Lcoh11,DES=="cens",select=c(TstopY,Tstopage,lex.id)) # select open episodes
   	  #  in next code: Lcoh12 replaced by Lcoh11
       text (Lcoh11$TstopY+1.0,
                Lcoh11$Tstopage,Lcoh11$lex.id,cex=0.7,adj=0)
   } else Lcoh12 <- "Lcoh12 not produced because more than 20 subjects were selected."
  legend ("topleft",legend=namstates,lty=1,col=colours,bg="white",cex=0.8)
  return (list(Lcoh11=Lcoh11,
               sub =subjectsID2,
               k = length(subjectsID2)))
}
