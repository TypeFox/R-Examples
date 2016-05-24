Biograph.Lexis <-
function (Bdata,Dlong) # ,Dlong,subjectsID,title)
{  # require (Epi)
  # Convert dates in years
   if (is.null(attr(Bdata,"format.date"))) {print ("Biograph.Lexis: format.date is missing (attribute of data)")}
  z<- check.par (Bdata)
  if (missing(Dlong)) 
       { print ("Getting data in long file format. Patience please.",quote=FALSE)
       	 Dlong <- Biograph.long(Bdata)
       	 Dlong2 <- Dlong$Depisode
       	 print ("Data in long format produced. Lexis continues.") } else Dlong2 <- Dlong # on input, Dlong is Dlong#Depisode
       	 
  if (!is.data.frame(Dlong2)) {stop ("Dlong$Depisode is not a data frame. Please check") }
  namstates <- attr(Bdata,"param")$namstates
   format.in <- attr(Bdata,"format.date")
   format.born <- attr(Bdata,"format.born")
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
  { print ("Biograph.Lexis.R: some durations are negative.")
    print (Dlong2[duration<0,])
    return
  }
 xx <- ifelse (as.character(Dlong2$DES)=="cens",as.character(Dlong2$OR),as.character(Dlong2$DES))
  Lcoh <- Lexis( id = Dlong2$ID,
               entry = list( per=en1 ),
               exit  = list( per=ex1, age=ex1-bt ),
               entry.status=as.character (Dlong2$OR),
               exit.status = xx,
               data=Dlong2,
               merge=TRUE)  
    rownames (Lcoh) <- 1:nrow(Lcoh)
    return (Lcoh)
}
