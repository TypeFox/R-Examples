TransitionAB <-
function (Bdata,transition,keep)
{ if (missing (transition)) stop ("TransitionAB: transition missing")
	f <- StateSpace (Bdata)
	namstates <- attr(Bdata,"param")$namstates
	if ((substring(transition,1,1)%in%namstates | substring(transition,1,1)=="*") ==FALSE) stop ("TransitionAB: at least one state not in state space. Try running StateSpace before this function")
	if ((substring(transition,2,2)%in%namstates | substring(transition,2,2)=="*") ==FALSE) stop ("TransitionAB: at least one state not in state space. Try running StateSpace before this function")
	if (missing(keep)) keep<-FALSE # do not remove the subjects that do not experience the transition
    format.in = attr(Bdata,"format.date")
    format.born <- attr(Bdata,"format.born")
	z<- check.par (Bdata)
	locpat <- locpath(Bdata)
	# Determine whether transition is from origin to destination state 
  # or from any origin to a destination state
  if (nchar(transition)!=2) 
   { print ("TransitionB: error in defining transition, use '*T' if you want destination T from any origin") 
     return
   }
  if (nchar(transition)==2 & substr(transition,1,1)=="*") ncasetrans <- 1 else ncasetrans <- 2  

# Determine, from Bdata$path, for each subject the position of the transition 
 pos <- sapply(Bdata$path,function(x) 
    {# if transition ="*N" and N is first state, remove first state
     if (transition==paste("*",substring(x,1,1),sep="")) 
            {aa=1; x <- substring(x,2,nchar(x))} else {aa=0}
     z1 <- grep(transition,x)
     pos <- ifelse (length(z1)==0,NA,
            ifelse (z1==1,nchar(unlist(strsplit(x,transition)[1]))+1,0))
            # Find all positions: which(strsplit(string, '')[[1]]=='a') (see pos.char.r)
     pos <- pos + aa
     if (substring(transition,2,2)=="*") pos <- ifelse (nchar(x)==pos,NA,pos) 
     return(pos)} )
 pos <- unname(pos)
 # if destination is *, there should be at least one state after the origin state
 
 # id <- ifelse (is.na(pos),NA,Bdata$ID) # identification number
 
  print (paste ("Age profile: Number of individuals with transition ",transition," is ",length(na.omit(pos)),sep=""),quote=FALSE)
  # = position of transition given by cmeanhar in Bdata$path
  if (keep==FALSE) 
     { Bdata <- subset(Bdata,!is.na(pos)) 
       pos <- na.omit(pos) } 
  # aa <- ifelse (ncasetrans==1,locpat+as.numeric(x[1]),locpat+as.numeric(x[1])+1)
  dates <- apply(cbind(pos,Bdata),1,function(x) 
             {z <- ifelse (ncasetrans==1,locpat+as.numeric(x[1]),locpat+as.numeric(x[1])+1)
              x[z]})
  dates <- as.numeric(dates)

  if (is.null(format.in)) stop("Biograph object: date format not specified. Check attributes of Biograph object. ")
  if (format.in=="age") {ages=dates
  	                     byear=as.numeric(Bdata$born)
  	                     years=rep(NA,nrow(Bdata))} else {
  y <- date_convert (dates,format.in=format.in,selectday=1,format.out="year",born=Bdata$born,format.born=format.born)
  years <- y
  y<- date_convert(as.numeric(Bdata$born),format.in=format.born,format.out="year",format.born=format.born)
  byear <- y
  ages <- years-byear  }
    
   cases <- c("From any origin to given destination","From given origin to given destination")
   cases <- paste ("Transition ",transition,": ",cases,sep="")
   return (list (case = cases[ncasetrans],
                 n = length(na.omit(ages)),
                 id=Bdata$ID,
                 pos = pos,
                 date = dates,
                 age = ages,
                 year=years,
                 cohort=byear))
}
