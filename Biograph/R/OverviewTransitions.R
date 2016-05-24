OverviewTransitions <-
function (Bdata,seq.ind,agetrans) {
# seq.ind produced by statesequence.ind.r
# agetrans produced by AgeTrans.r
 # trans[ix,io,id] = Transitions by age [ix], origin [io] and destination [id]
#   origin = state before j-th transition: seq.ind[i,j]
#   destin = state after j-th transition: seq.ind[i,j+1]
# trans[i,j] is from seq.ind[i,j] to seq.ind[i,j+1]     
#    (eg Bdata[6,] ; ages[6,]; events[6,] ; seq.ind[6,]
# input
#   -   tstate   Global variable
# Direct transitions and censored cases
   z<- check.par(Bdata)
   namstates <- attr(Bdata,"param")$namstates
  numstates <- length (namstates)
  iagelow <- attr(Bdata,"param")$iagelow
  iagehigh <- attr(Bdata,"param")$iagehigh
  nage <- attr(Bdata,"param")$nage
  format.born <- attr(Bdata,"format.born")
   if (missing(seq.ind)) 
       { print ("Calling function seq.ind",quote=FALSE)
       	seq.ind <- Sequences.ind (Bdata$path,namstates) }
   if (missing(agetrans)) 
      {  print ("Calling function AgetTrans",quote=FALSE)
      	 agetrans <- AgeTrans (Bdata)  }
  ages <- agetrans$ages
  agecens <- agetrans$agecens
  nsample <- nrow(Bdata)
  agelist2 <- c(iagelow:iagehigh) 
  zz <- length(agelist2)  
  ns <- nchar(Bdata$path)           
  trans <- array(0,c(zz,numstates,numstates+1))
  zz1 <- zz - 1                 
  dimnames(trans) <- list (Age = c(0:zz1),origin=namstates,destination=c(namstates,"censored"))
  if (!is.null(attr(Bdata,"format.date"))) 
   {format.in <- attr(Bdata,"format.date")} else  {z <- check.par (Bdata)}
  if (attr(Bdata,"format.date") != "days" & max(agecens) > 150) warning("OverviewTransitions: age at censoring exceeds 150. Please check.")
   
  for (i in 1:nsample) {
  	yb <- date_convert (Bdata$born[i],format.in=format.born,format.out="year",born=Bdata$born[i],format.born=format.born)
    y <- date_convert (d=Bdata$end[i],format.in=format.in,format.out="year",born=Bdata$born[i],format.born=format.born)
    agecens <- ifelse (ns[i]==1,trunc(y-yb)-iagelow+1,ifelse(y==ages[i,(ns[i]-1)],NA,trunc(y-yb)-iagelow+1))
    
    if (ns[i] > 1) 
      { for (j in 1:ns[i]-1) 
       {  age_at_trans <- trunc(ages[i,j])-iagelow+1   # + 1 because of vector definition
#         occupx[i,ix] <<-  seq.ind         
         trans[age_at_trans,seq.ind[i,j],seq.ind[i,j+1]] <-  trans[age_at_trans,seq.ind[i,j],seq.ind[i,j+1]] + 1
       }}
    if (!is.na(agecens)) trans[agecens,seq.ind[i,ns[i]],numstates+1] <- trans[agecens,seq.ind[i,ns[i]],numstates+1] + 1
  }                                                                                     
 #occupxt <- tstate
  
  # Total number of direct transitions and censored cases 
  n1 <- numstates + 1 
  n2 <- numstates + 2
  n3 <- numstates + 3
  Ttrans <-array(0,c(numstates+1,numstates+3))
  Ttrans[1:numstates,1:n1] <- apply(trans,c(2,3),sum)
  Ttrans[,n2] <- Ttrans[,n1]  # move censored cases to last column
  Ttrans[,n1] <- 0
  Ttrans[,n1] <- rowSums(Ttrans[,-n2])
  Ttrans[n1,] <- colSums(Ttrans)
  Ttrans[,n3] <- Ttrans[,n2] +Ttrans[,n1]
  dimnames(Ttrans) <- list(Origin=c(namstates,"Total"),Destination=c(namstates,"Total","Censored","TOTAL"))
  #trans_possible <<- trans_possible
 # Mean age at transition 
  sweep(trans,c(2,3),colSums(trans),"+")[1,,] #total number of transitions by origin and destination
  agedistrib <- sweep(trans,c(2,3),colSums(trans),"/")  # age distribution
  aged <- array(0,c(iagehigh-iagelow,numstates,numstates+1))
  for (ix in 1:(nage-1)) aged[ix,,] <- agedistrib[ix,,] * (agelist2[ix]-0.5) #  was iagelow:(iagehigh-1)  7 Aug 08
  meanage <- apply(aged,c(2,3),sum)
  dimnames(meanage) <- list(Origin=namstates,Destination=c(namstates,"censored"))
  meanage <- round(meanage,2)  
       
 return (list ( Ttrans = Ttrans,
                meanage = meanage))
 }
