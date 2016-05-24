SamplePath <-
function (Bdata,subjectsID)
# Produces sample path for selection of subjects
#       The record numbers of the selected subjects are in vector "subjects"
#       e.g. subjects <- c(1,40,765,5320)
# NEEDS Parameters (first run Parameters)
# Output: prints sample paths to text file; no other output
{  namstates <- attr(Bdata,"param")$namstates
	if (missing(subjectsID)) 
	  { print ("SamplePath: Path for random sample of 5 subjects",quote=FALSE)
	  	subjectsID <- sample(Bdata$ID,5,replace=FALSE) }
  z <- check.par (Bdata)
  locpat <- attr(Bdata,"param")$locpat
  seq.ind <- Sequences.ind (Bdata$path,namstates)
 ages <- AgeTrans(Bdata)$ages  
 calendar <- date_b (Bdata,selectday=1,format.out="%d%b%y")
# if format.in = age and Bdata$born not provided, skip 

  nn <- length(subjectsID)
  # Remove IDs that are not in sample (in survey data)
  # Check whether subjectsID are in ID
   z<- subjectsID %in% Bdata$ID
   if (TRUE %in% z==FALSE) 
      {print ("print.Samplepath: wrong subjectsID. Check ID numbers.")
       subjectsID <- Bdata$ID[1]
      } 
  for (k in 1:nn) 
  { i <- which (Bdata$ID == subjectsID[k])
    if (length(i) == 0) subjectsID[k] <- NA
  }
  subjectsID <- subset(subjectsID,!is.na(subjectsID))
  nn <- length(subjectsID)
  ID <- vector (mode="numeric",length=nn)
  born <- vector (mode="character",length=nn)
  path <- vector ("list",nn)
  info <- vector ("list",nn)
  zc <- c(1:nrow(Bdata))
  for (k in 1:nn) 
  { # i<- as.numeric(rownames(survey)[survey$ID==subjectsID[k]]) Replaced 7/5/2010
     i <- zc[Bdata$ID==subjectsID[k]]
       # subjectsID = ID; get associated line number
     ns <- nchar(Bdata$path[i]) 
    nss <- ns
    nss1 <- nss + 1 # = number of transitions + 2 
    nss2 <- nss-1 
    cmc_trans <- vector(mode="numeric",length=nss1)
    cal_trans <-  vector(mode="numeric",length=nss1)
    cmc_dur <- vector(mode="numeric",length=nss)
    n78 <- locpat+nss # gives cmc at last transition before censoring  
    cmc_trans[1] <- Bdata$start[i]
    cal_trans[1] <- calendar$start[i] 
    if (nss == 1) 
      { cmc_trans[2] <- Bdata$end[i]
      	cal_trans[2] <- calendar$end[i] } else 
      { cmc_trans[2:nss1] <- c(Bdata[i,(locpat+1):(n78-1)],Bdata$end[i])
      	cal_trans[2:nss1] <- c(calendar[i,(locpat+1):(n78-1)],calendar$end[i]) }
    cmc_trans <- unlist (unname(cmc_trans))
    cal_trans <- unlist (unname(cal_trans))
    for (j in 1:nss) cmc_dur[j] <- (cmc_trans[j+1]-cmc_trans[j])
    cmc_dur[nss1] <- NA
      if (!is.null(attr(Bdata,"format.date"))) 
   {format.in <- attr(Bdata,"format.date")} else  stop('p.SamplePath: format.date not specified (see attr(Bdata,"format.date"))')
   format.born <- attr(Bdata,"format.born")
    y <- date_convert (Bdata$born[i],format.in=format.in,format.out="year",format.born=format.born) 
    born2 <- y 
    y <- date_convert (Bdata$start[i],format.in=format.in,format.out="year",format.born=format.born)
    age_start <- y-born2
    y <- date_convert (Bdata$end[i],format.in=format.in,format.out="year",format.born=format.born)
    age_cens <- y-born2
 #   ntimeunit <- ifelse (attr(Bdata,"format.date")=="month",12,ifelse (attr(Bdata,"format.date")=="year",1,1))
 #   age_start <- (Bdata$start[i]-Bdata$born[i])/ntimeunit
 #   age_cens <- (Bdata$end[i]-Bdata$born[i])/ntimeunit
    ID[k] <- Bdata$ID[i]
    born[k] <-  paste("Subject ID = ",ID[k], "  Date of birth ",round(Bdata$born[i],2)," (",calendar$born[i],")",sep="")
    
  #  print (born[k],quote=FALSE)
  #  print ("Observation window:",quote=FALSE)
  #  print (paste("Start ",round(Bdata$start[i],2),"  (",calendar$start[i],")",sep=""),quote=FALSE)
  #  print (paste("End   ",round(Bdata$end[i],2)  ,"  (",calendar$end[i],")",sep=""),quote=FALSE)
    dur <- paste("Total duration of observation ",round(Bdata$end[i]-Bdata$start[i],2),
                 " ",attr(Bdata,"format.date"),"s",sep="")
 #   if (nss ==1)  <- c(round(age_start,2),round(age_cens,2)) else 
          age55 <- c(round(age_start,2),c(round(ages[i,1:nss2],2),round(age_cens,2)))
    
    path[[k]] <- data.frame(Episode=1:nss1,State=c(namstates[seq.ind[i,1:nss]],"Cens"),EntryDate1=c(round(cmc_trans,2)),EntryDate2=cal_trans, EntryAge=age55,Durat=round(cmc_dur[1:nss1],2),OR=c(0,seq.ind[i,1:nss]),DE=c(seq.ind[i,1:nss],0))
   # print (path[[k]])
    rownames(path)<- NULL
    TAB_samplepath <-noquote(format(path[[k]],justify="right"))
  #  print (TAB_samplepath)

    info[[k]] <- list(ID=ID[k],born=born[k],path=path[[k]])
  }

return(info)
}
