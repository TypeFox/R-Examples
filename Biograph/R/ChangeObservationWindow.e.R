ChangeObservationWindow.e <-
function (Bdata,entrystate,exitstate)
{ #  Check whether Parameters was called
  z<- check.par(Bdata)
  	namstates <- attr(Bdata,"param")$namstates
  # Check presence of entrystate and exitstate in Bdata. Remove records with missing entry- or exitstate
  # This section copied from function TransitionAB()
   if (entrystate%in%namstates & (exitstate%in%namstates|is.na(exitstate))) z=3 else stop("ChangeObservationWindow.e : Entry- or exitstate not part of state space")
   if (!is.na(entrystate))
   {  transition <- entrystate
      pos <- sapply(Bdata$path, function(x) {
         z1 <- grep(transition, x)
         pos <- ifelse(length(z1) == 0, NA, ifelse(z1 == 1, nchar(unlist(strsplit(x, 
              transition)[1])) + 1, 0))
       })
     pos <- unname(pos)
     Bdata6 <- subset (Bdata,!is.na(pos))
     attr(Bdata6,"format.date") <- attr(Bdata,"format.date")
     attr(Bdata6,"format.born") <- attr(Bdata,"format.born")
     attr(Bdata6,"param") <- attr(Bdata,"param")
     Bdata <- Bdata6
    }
  print (paste(nrow(Bdata)," observations with entrystate ",entrystate,sep=""),quote=FALSE)
  if (!is.na(exitstate))
  { transition <- exitstate
    pos <- sapply(Bdata$path, function(x) {
      z1 <- grep(transition, x)
      pos <- ifelse(length(z1) == 0, NA, ifelse(z1 == 1, nchar(unlist(strsplit(x, 
            transition)[1])) + 1, 0))
      })
     pos <- unname(pos)
     Bdata3 <- subset (Bdata,!is.na(pos))
     attr(Bdata3,"format.date") <- attr(Bdata,"format.date")
     attr(Bdata3,"format.born") <- attr(Bdata,"format.born")
     attr(Bdata3,"param") <- attr(Bdata,"param")
     Bdata <- Bdata3
   }  
  print (paste(nrow(Bdata)," observations with entrystate ",entrystate," and  exitstate ",exitstate,sep=""),quote=FALSE)
  locpat <- locpath(Bdata)
  maxtrans <- (ncol(Bdata)-locpat) 
  windowe <- function (data,refstate)
  { entry88 <- function (x,state)
    {   # Get for a state following another state the position in path and the date of occurrence
        #  x[locpath(Bdata)] is the state sequence path
        #  x[length(x)] is the first state occupied in current sequence (=low)
        #  x[locpath(Bdata)-1] is the last state occupied (= ns) REMOVED
        
    	lowx <- as.numeric(x[length(x)]) # see below
    	if (is.na(lowx))
    	 { return (list(loctrans = NA,
    	 	            date = NA))} else
    	# Is state in substring that starts at position 'lowx' [= x[length(x)]] in path?   
    	{ yn <- pos.char(x[locpat],state)
    	  if (is.na(yn))  # state in not in path
    	    {return (list(loctrans = NA,
    	    	            date = NA))} else
         { # State is in substring of path
           # Get the substring
             ns <- nchar(x[locpat])
             kk <- substr(x[locpat],x[length(x)],ns) # x[(locpat-1)])
           # Determine position of state (loc) in substring and date at entry into state
             loc <- grep(state,stringf(kk))[1]  # loc=first location of (reference) state
          # z98 <- ifelse (is.na(loc),ns-as.numeric(x[length(x)])+1,(z99 -1 + as.numeric(x[length(x)])) )  # as.numeric(x[(locpat-1)])
          # Position of state in path
           locp <- ifelse (loc==1,lowx,lowx + loc -1)
           date <- ifelse (loc==1,x[3],x[(locpat+locp-1)]) # date at entry in state
            # if loc==1, date is Bdata$start
           # if (is.na(state)) date <- x[(locpat+z98)]
           return(list(loctrans = locp,    #  location of transition (entry or exit)
                       date=date))   #  date of transition
         } 
       }
    } # end entry88
    
    if (is.na(refstate)) # user gives NA as reference state
     { loc <- nchar(data[,locpat])   #    data[,(locpat-1)]
       date <- data[,4]
     }  else
     { date90 <- apply (data,1,function(x) entry88(x,refstate))
       xx<- sapply(date90,function(x) {x[1]})
       loc <- unname(unlist(xx)) #  location of state at entry in observation
       xx<- sapply(date90,function(x) x[2])
       date <- as.numeric(unname(unlist(xx)) )# date at entry in observation 
     }
    return (list (loc=loc,
                  date=date))
  }  # end windowe
  print ("Creating new Biograph object. Patience please . . . ")
  low1 <- rep(1,nrow(Bdata))
  entry <- windowe (cbind(Bdata,low1),entrystate)   # GLHS low[1] = 2
  exit <- windowe (cbind(Bdata,low1=entry$loc),exitstate)
 # exit$loc <- exit$loc + entry$loc - rep(1,length(exit$loc))  # CHANGED 10_11_2010
   # nn <- nrow(Bdata[!is.na(exit$loc),]) 
   
# Create a new Biogaph data file: subset with observations in the interval
  entryloc <- entry$loc
  exitloc <-  exit$loc
  #cmc_select <- subset(cmc,!is.na(entrydate))
  entrydate <- entry$date  # STARTING DATE for subjects with entry event
  exitdate <-  ifelse (is.na(exit$date),Bdata$end,exit$date)  # idem
  #  if exit (censoring) can be at Bdata date
  Bdata2 <- Bdata
  for (i in 1:nrow(Bdata))
   { Bdata2$path[i] <- ifelse (is.na(entry$loc[i])," ",ifelse (is.na(exit$date[i]),entrystate,substr(Bdata$path[i],entry$loc[i],exit$loc[i])))}  
  print ("Continues creating new Biograph object. Patience please . . . ")
  for (i in 1:nrow(Bdata))
  {  ns <-  nchar(Bdata2$path[i])  
     Bdata2$start[i] <- ifelse (entry$loc[i]==1,Bdata$start[i],Bdata[i,(locpat+entry$loc[i]-1)])
     nn <- ifelse (ns==1,1,ns-1)
     if (nn==1) {Bdata2[i,(locpat+1)] <- Bdata[i,(locpat+entry$loc[i])]} else
           {Bdata2[i,(locpat+1):(locpat+nn)] <- Bdata[i,(locpat+entry$loc[i]):(locpat+exit$loc[i]-1)]}
     #   ifelse (is.na(exit$loc[i]),NA,Bdata[i,(locpat+entry$loc[i]):(locpat+exit$loc[i]-1)])
     Bdata2[i,(locpat+ns):(locpat+maxtrans)] <- NA  
     Bdata2$end <- exitdate
   }
   Bdata2 <- Bdata2[Bdata2$path!=" ",]
   attr(Bdata2,"format.date") <- attr(Bdata,"format.date")
   attr(Bdata2,"format.born") <- attr(Bdata,"format.born")
   param <- Parameters(Bdata2)
   attr(Bdata2, "param") <- param
   print("A Biograph object with new observation window is returned.",quote = FALSE)
  return (Bdata2)
 }
