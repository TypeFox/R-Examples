Biograph.long <-
function (Bdata)
{ z <- check.par (Bdata)
  # Remove records with end < start
  Bdata2 <- subset (Bdata[Bdata$end-Bdata$start >= 0,]) # was: Bdata
  attr(Bdata2,"format.date") <- format.in <- attr(Bdata,"format.date")
  if (is.null(attr(Bdata,"param"))) print ("Biograph.long: Parameters missing. Run Parameters . . . . ",quote=FALSE)
  tmat <- attr(Bdata,"param")$tmat
 # attr(Bdata2,"trans") <- attr(Bdata,"trans")
 # if (!is.null(attr(Bdata2,"trans"))) tmat <- attr(Bdata2,"trans") else 
 #           #  tmat <- ex_Bdata$tmat 
 #    {print ("Biograph.long: Parameters missing. Biograph runs Parameters . . . . ",quote=FALSE)
 #     ex_Bdata <- Parameters(Bdata2)
 #     attr(Bdata2,"trans") <- ex_Bdata$tmat
 #     tmat <- ex_Bdata$tmat
 #    }
  print (". . . . .  Creating long format  . . . . . .",quote=FALSE)
   locpat <- locpath(Bdata2)
  print (". . . . .  running reshape  . . . . . ",quote=FALSE)
  # require (reshape)
  nn99 <- locpat + max(nchar(Bdata2$path))-1
  zx <- reshape (Bdata2,idvar="ID",varying=list(c(3,(locpat+1):nn99,4)),
     v.names="date",direction="long",drop=NULL)
  print (" . . . . Sort data in long format  . . . . ",quote=FALSE)
  zx2 <- zx[!is.na(zx$date),]
  D <- zx2[do.call("order",list(zx2$ID,zx2$date)),] # sort by 2 variables
  print (" . . . . Adjust long format for survival package etc . . .  ",quote=FALSE)
  D$OD <- substr(paste("B",D$path,"Ce",sep=""),D$time,D$time+1)
  D$ns <- nchar(D$path)
  D$time <- ifelse (D$time==max(D$time),D$ns+1,D$time)   # time = line number of episode in trajectory (first= from birth; last = open to censored)
  D$OD[D$time==D$ns+1]<- "cens"
  D$Tstart <- rep(0,nrow(D))
  D$Tstart[2:nrow(D)] <- ifelse (as.numeric(D$date[2:nrow(D)]) ==as.numeric(D$born[2:nrow(D)]),D$born[2:nrow(D)],D$date[2:nrow(D)] -diff(D$date,lag=1))
  if (class(D$date)=="Date") D$Tstart <- as.Date(D$Tstart,origin="1970-01-01")
  D$Tstop <- D$date

# ----------------------------------
  ns <- nchar(Bdata2$path)
  max_ns <- max(ns)
 # D8 <- subset(D,D$Tstart != D$Tstop & D$Tstart >= D$born)  # remove first episode (birth to entry in first state)  was: D$date!=D$born & 
  D8<- subset (D,D$time<=(max(D$ns)+1))  # & D$OD!="")
  D8 <-subset (D8,D8$time!=1)
  D8 <- subset (D8,!(D8$date==D8$Tstart & D8$time==(max_ns+1)))
   #  in Jan 2012, D8$born was replaced by D8$Tstart (running format.dat=age)
   # and ChangeObservationWindow.t to start at age 10 (GLHS)
  D_negative <- subset (D8,D8$Tstop-D8$Tstart < 0) # = empty
  D <- D8
    # Check whether duration is nonnegative (for Lexislines): select records with neg values
    # result should be same as sum(Bdata$ns)
  if (nrow(D)!=sum(ns)) warning (paste("Biograph.long: Number of records in long format differs from total number of episodes in Biograph object: ",nrow(D)," versus ",sum(ns),sep=""))
  # At birth: state is B or first state of state sequence
  CC <- "B"
  CC <- substr(D$path,1,1)
  D$OR <- ifelse (D$Tstart==0,CC,ifelse (D$time >D$ns, substr(D$path,D$ns,D$ns), substr(D$path,(D$time-1),(D$time-1))))  
  D$DES <-   ifelse (D$time > D$ns, "cens",substr(D$path,(D$time),(D$time)))

  # D$DES <- ifelse (D$time > D$ns,"cens",substr(D$path,(D$time),(D$time)))
  D$status <- ifelse (D$DES=="cens",0,1)
  zloc=9
  ncovariates <- attr(Bdata,"param")$ncovariates
   namstates <-   attr(Bdata,"param")$namstates
  D$trans <- apply(D,1,function (x) {ifelse (x[ncovariates+zloc+2]=="cens", 
           grep(x[ncovariates+zloc+1],namstates),
           tmat[grep(x[ncovariates+zloc+1],namstates),grep(x[ncovariates+zloc+2],namstates)])})
# D$trans changed october 2011
  code88 <- ifelse(format.in=="CMC",12,1)
  D$birthdate <- D$born
  if (class(D$born)=="Date") D$birthyear <- Date_as_year (D$born) else 
     {{ if (attr(Bdata,"format.date")=="CMC"|attr(Bdata,"format.date")=="cmc") D$birthyear <- date_convert (D$born,format.in="CMC",format.out="year",format.born=attr(Bdata,"format.born")) else
     	  { if (attr(Bdata,"format.date")=="year"|attr(Bdata,"format.date")=="age"|attr(Bdata,"format.date")=="numeric")  D$birthyear <- D$born else D$birthyear <- D$born}
     }}	  
  if (attr(Bdata,"format.date") == "age") D$born <- rep(0,nrow(D)) 
  D$Tstarta <- (D$Tstart-D$born)/code88
  D$Tstopa  <- (D$Tstop-D$born)/code88
  
  De <- subset(D,D$time!=1)  # remove transition from birth to first state
  De$time <- De$time-1
  Depisode = data.frame (ID=De$ID,OR=De$OR,DES=De$DES,Tstart=De$Tstart,
                Tstop=De$Tstop,status=De$status,trans=De$trans,
                birthdate=De$birthdate,birthyear=D$birthyear,De[,(3:(2+ncovariates))],
                born=De$born,OD=De$OD,Episode=De$time,Tstarta=D$Tstarta,Tstopa=D$Tstopa)
                         
  attr(D, "param") <- attr(Bdata,"param")  
  attr(Depisode, "param") <- attr(Bdata,"param") 
  attr(D, "format.date") <- attr(Bdata,"format.date")   
  attr(Depisode, "format.date") <- attr(Bdata,"format.date") 
  attr(D, "format.born") <- attr(Bdata,"format.born")   
  attr(Depisode, "format.born") <- attr(Bdata,"format.born") 
  
  return (list (Devent = D,
                Depisode = Depisode))
}
