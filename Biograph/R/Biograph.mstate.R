Biograph.mstate <-
function (Bdata)
# Prepares data format for mstate (Putter)
{ 	z <- check.par(Bdata)
	namstates <-  attr(Bdata,"param")$namstates
  numstates <- length(namstates)
   tmat <- attr(Bdata,"param")$tmat
   if (is.null(attr(Bdata,"param"))) print ("Biograph.mstate: Parameters missing. Run Parameters . . . . ",quote=FALSE)
 covnames <- colnames(Bdata)[5:(locpath(Bdata)-3)]
 # 1. Remove intrastate transitions  ========
  z<- !is.na(diag(attr(Bdata,"param")$tmat))
  if (TRUE %in%z) # at least one diagonal element is not missing
  {  print ("     Biograph.mstate: calls function Remove_intrastate  . . . ")
  dd <- Remove.intrastate (Bdata)
  Bdata2 <- dd
  param <- Parameters(Bdata2)
  tmat <- param$tmat
  } else {Bdata2 <- Bdata 
          attr(Bdata2,"param") <- attr(Bdata,"param")  }
 # attr (Bdata2,"trans") <- tmat
  # 2. Produce long format  ===================
  z<- Biograph.long (Bdata2)
  Dlong  <- z$Depisode                                                        

  # 3. Create mstate format from long format ==
  # namstates <- unique(Dlong$OR)
  # numstates <- length(namstates)
  # tmat <- attr(Dlong,"trans")
  inamstates <- vector (mode="numeric",length=numstates)
  for (idd in 1:numstates)
  { inamstates[idd] <- grep(namstates[idd],namstates)  
  }
  from <- apply(Dlong,1, function(x) grep(x[2],namstates))  
  to <- apply (Dlong,1, function(x) ifelse (x[3]=="cens","cens",grep(x[3],namstates)))
  Dlong$from <- from
  Dlong$to <- to

# Create mstate format (D2) from long format (Dlong)
nrows = nrow(Dlong)
D2 <- Dlong[1:nrows,]

# in tmat: number of destinations by origin
# tmat: diagonal set from NA to zero (for mvna). Replace
diag(tmat) <- NA  # DIAGONAL OF TMAT MUST BE NA, NOT 0 !!!!!!!! (DEC 2010)
nt <- apply(tmat,1,function(x) length(na.omit(x))) 
# Obtain for each transition, origin and destination (ie position on trans-matrix)
or <- des <- vector (mode="numeric",length=max(tmat,na.rm=TRUE))
for (i in 1:max(tmat,na.rm=TRUE)) # transition numbers (values of matrix elements)
 {zk <- which (tmat==i)
  or[i] <- zk - nrow(tmat) * trunc((zk-1)/nrow(tmat))
  des[i] <- trunc((zk-1)/nrow(tmat))+1  
 }
# number of destinations by record in Dlong
nidd <- nt[Dlong$from[1:nrows]] 
# Number of records in D2 data frame: tt (records in Dlong: nrow(Dlong))
tt <- sum(nt[Dlong$from[1:nrows]])
# Get indicators in mstate long format (tt records)
zx <- array(0,dim=c(tt,4))
zx[,1] <- 1:tt
jk <- 0
for (i in 1:nrows) 
 {if (nidd[i] > 0 )
 	{   zx[(jk+1):(jk+nidd[i]),2]=i
 		zx[(jk+1):(jk+nidd[i]),3]<- Dlong$from[i]
 		zx[(jk+1):(jk+nidd[i]),4]<- as.numeric(ifelse(Dlong$to[i]=="cens",0,Dlong$to[i])) # includes "cens"
 		jk<-jk+nidd[i]
  } }

# Get for each episode (transition), the possible destinations: trans.possible
# in diagonal = 0, trans.possible = 0, which creates problem later
trans.possible <-array(NA,dim=c(nrows,numstates)) 
for (i in 1:nrows)
  { if (nidd[i] >0 ) ib <- as.numeric(na.omit(tmat[Dlong$from[i],])) else ib <- NA
  	trans.possible[i,1:length(ib)] <- ib
  }
trans.possible[trans.possible==0] <- NA
 
#  Replace destination with other possible destination
D2[zx[,1],]<-Dlong[zx[,2],]
z <- t(trans.possible) # if z is zero, D2$trans = 0 and des=missing   
     #  hence length of des[D2$trans] zx[,4]  differ 
     #  length(z) = nrow(Dlong) * numstates
D2$trans <- z[!is.na(z)] 
#$ status = 1 if destination derived from trans is "to"
D2$status <- 0
D2$status[D2$trans>0] <- ifelse (des[D2$trans]==zx[,4][D2$trans>0],1,0) # event or no event
D2$status <- ifelse (des[D2$trans]==zx[,4][D2$trans>0],1,0)
D2$DES <-namstates[des[D2$trans]] 
D2$to <-inamstates[des[D2$trans]] 

attr(D2,"trans") <- attr(Bdata2,"param")$tmat
attr(D2, "param") <- attr(Bdata2,"param")                                                         
class(D2) <- c("msdata", "data.frame")
return (D2)
}
