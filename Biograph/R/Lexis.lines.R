Lexis.lines <-
function (Bdata,Dlong,subjectsID,title)
  # NOTE: Dlong = Dlong$Depisode
  # Bdata =GLHS.yr
 {   # Convert dates in years
  if (missing(subjectsID)) subjectsID <- sample(Bdata$ID,5,replace=FALSE)
  if (missing(title)) title <- "Title missing"
  namstates <- attr(Bdata,"param")$namstates
  z<- check.par (Bdata)
  format.in <- attr(Bdata,"param")$format.date
  format.born <- attr(Bdata,"param")$format.born
  # attr(Bdata,"format.date") should be 'year'
  # if date or birth is CMC, than change it in years
  # Convert the dates into years (also date of birth)
  if (format.in!="year")  Bdata <- date_b (Bdata,selectday=1,format.out="year")
   format.in <- attr(Bdata,"param")$format.date
  format.born <- attr(Bdata,"param")$format.born
   if (missing(Dlong)| attr(Dlong,"format.date")!="year") 
       # get long format of Bdata (to get Dlong with same date format as Bdata)
       { print ("Lexis.lines: Getting data in long file format. Patience please.",quote=FALSE)
       	 D <- Biograph.long(Bdata)
       	 Dlong2 <- D$Depisode
       	 attr(Dlong2,"format.date") <- attr(Bdata,"format.date")
       	 attr(Dlong2,"param") <- attr(Bdata,"param")
       	 print ("Data in long format produced. Lexis continues.") } else Dlong2 <- Dlong # on input, Dlong is Dlong#Depisode
   if (attr(Bdata,"format.date")!=attr(Dlong2,"format.date")) stop (paste("Date format in Dlong is not year. It is: ",attr(Dlong2,"format.date"),sep="") )  	 
  if (!is.data.frame(Dlong2)) {stop ("Dlong$Depisode is not a data frame. Please check") }
  
  if (!is.null(attr(Bdata,"format.date"))) 
   {format.in <- attr(Bdata,"format.date")} else  {print ("Lexisines.episodes: format.date is missing (attribute of data)")}
   format.born <- attr(Bdata,"format.born")
   y <- date_convert (Dlong2$Tstart,format.in=format.in,format.out="year",format.born=format.born)
   Dlong2$TstartY <- y
   y <- date_convert (Dlong2$Tstop,format.in=format.in,format.out="year",format.born=format.born)
   Dlong2$TstopY <- y
   y <-date_convert (Dlong2$born,format.in=format.born,format.out="year",format.born=format.born) 
   bt <- y
   Dlong2$Tstartage <- Dlong2$TstartY - Dlong2$born
   Dlong2$Tstopage <- Dlong2$TstopY - Dlong2$born
         
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
age1 <- en1-bt
age2 <- ex1-bt
# Create input date file
nsample <- nrow(Bdata)

print ("Producing long format for ggplot2. Patience please.")
selectionID <- which(Bdata$ID%in%subjectsID)
jj = 0
for (i in selectionID)    #  1:nsample)
  { jj <- jj + 1
  	zd50 <- 100*trunc(jj/100)
  	if (jj==zd50) print (paste(jj," records of ",length(selectionID)," processed",sep=""),quote=FALSE)
   z1 <- Dlong2[Dlong2$ID==Bdata$ID[i],]
   z2 <- Dlong2[Dlong2$ID==Bdata$ID[i]&Dlong2$Tstop==Bdata$end[i],]
   z2$Tstart <- z2$Tstop
   z2$TstartY <- z2$TstopY
   z2$Tstarta <- z2$Tstopa

  	dat4 <- rbind (z1,z2)
  	dat4$locIDx <-  date_convert (Bdata$end[i],format.in=format.in,format.out="year",format.born=format.born) 
  	dat4$locIDy <- date_convert (Bdata$end[i],format.in=format.in,format.out="age",born=Bdata$born[i],format.born=format.born) 
  	#(Bdata$end[i]-Bdata$born[i])
  	if (jj == 1)
  	  { data <- dat4} else
  	  { data <- rbind (data,dat4)}
  }
# data for plotting: "data" 
data$state <- factor (data$OR,levels=namstates,labels=namstates)

data2 <- subset(data,data$ID%in%subjectsID)
 nyear <- 5 
  if (max(na.omit(ex1-bt))-min(na.omit(en1-bt))<20) nyear <- 1  
  AgeLow <- nyear*trunc(min(na.omit(en1-bt))/nyear)
  AgeHigh <- nyear*trunc(max(na.omit(ex1-bt))/nyear+1)
  PerLow <- nyear* trunc(min(na.omit(en1)/nyear))
  PerHigh <- nyear* (trunc(max(na.omit(ex1)/nyear))+1)
  PerHigh[PerHigh-PerLow < AgeHigh-AgeLow] <- PerLow + AgeHigh - AgeLow
  AgeHigh[AgeHigh-AgeLow < PerHigh-PerLow] <- AgeLow + PerHigh - PerLow

colours <- c("red","green","orange","purple","brown","orange")
colours <- c("green","orange","red","yellow","orange","blue")
# -----------------------------------------------------------
# to avoid: no visible binding for global variable
# http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
TstartY <- Tstarta <- state <- ID <- locIDx <- locIDy <- NULL 
# ----------------------------------------------------------
lex3 <- ggplot(data2,aes(x=TstartY,y=Tstarta,colour=state))  # colour p. 48
#lex3 <- ggplot(data2,aes(x=x,y=y,colour=state))
p<- lex3+geom_line(aes(group = ID),size=1) + coord_equal()+scale_x_continuous(breaks=seq(PerLow,PerHigh+1,by=nyear))+scale_y_continuous(breaks=seq(AgeLow,AgeHigh+1,by=nyear)) # p. 108   title:p. 143
#p = p+geom_line (aes(x=1980),colour="yellow",size=1.5)
# print ID at end of lifeline
p2 <- p +geom_text(aes(label=ID,x=locIDx,y=locIDy,hjust=0, vjust=0),size=2.5,colour="yellow")
p2 <- p2 + xlab("Calendar year")+ylab("Age")
p3<- p2+ ggtitle(title)
p4 <- p3+theme(plot.title=element_text(size=11))
p5 <- p4+scale_colour_manual(values=colours)+theme(plot.background=element_rect(fill="lightskyblue1",colour=NA),
  panel.background=element_rect("black"),
  axis.text.x=element_text(colour="black"),
  axis.text.y=element_text(colour="black"),
  axis.title.x=element_text(colour="darkgreen",face="bold"),
  axis.title.y=element_text(colour="darkgreen",face="bold",angle=90))
print (p5)
return (list (subjectsID=subjectsID,
             p=p5))
}
