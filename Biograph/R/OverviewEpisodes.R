OverviewEpisodes <-
function (Bdata,seq.ind) {
  z <- check.par (Bdata)
  namstates <- attr(Bdata,"param")$namstates
  numstates <- length (namstates)
  nsample <- nrow(Bdata)
  locpat <- locpath(Bdata)
  ns <- nchar (Bdata$path)
  if (missing(seq.ind)) 
       { print ("Calling function seq.ind",quote=FALSE)
       	seq.ind <- Sequences.ind (Bdata$path,namstates) }

    # Open and closed episodes
  # no transition: determine state at onset of observation
  NumEpisodes <- array(0,c(numstates+1,5))
  for (i in 1:nsample) {
    nns <- nchar(Bdata$path[i])
    nns2 <- nns-1
  #     LRO and LO
    if(nns==1) NumEpisodes[seq.ind[i,1],1] <- NumEpisodes[seq.ind[i,1],1] + 1 
    else {NumEpisodes[seq.ind[i,1],2] <- NumEpisodes[seq.ind[i,1],2] + 1 
  # RO
          NumEpisodes[seq.ind[i,nns],3] <-  NumEpisodes[seq.ind[i,nns],3] + 1
  # Closed
          if (nns>2) for(j in 2:nns2) NumEpisodes[seq.ind[i,j],4] <- NumEpisodes[seq.ind[i,j],4]+1
          }
   }
  NumEpisodes[,5] <- rowSums(NumEpisodes)
  NumEpisodes[numstates+1,] <- colSums(NumEpisodes)
  dimnames(NumEpisodes) <- list(Episode=c(namstates,"Total"),Type=c("LROpen","LOpen","ROpen","Closed","Total"))
  
    DurEpisodes <- array(0,c(numstates+1,5))  
  for (i in 1:nsample) {
    nns <- nchar(Bdata$path[i])
    nns2 <- nns-1
    nns77 <- locpat+nns-1    # number of last transition before censoring
    if(nns==1) {DurEpisodes[seq.ind[i,1],1] <- DurEpisodes[seq.ind[i,1],1] + Bdata$end[i]-Bdata$start[i]} else 
       {DurEpisodes[seq.ind[i,1],2] <- DurEpisodes[seq.ind[i,1],2] + Bdata[i,(locpat+1)]-Bdata$start[i] 
          DurEpisodes[seq.ind[i,nns],3] <-  DurEpisodes[seq.ind[i,nns],3] + Bdata$end[i]-Bdata[i,nns77]
          if (nns>2) for(j in 2:nns2) { 
              jj21 <- locpat + j -1
              jj22 <- jj21 + 1 
              DurEpisodes[seq.ind[i,j],4] <- DurEpisodes[seq.ind[i,j],4]+ Bdata[i,jj22] - Bdata[i,jj21] }
          }
   }    
  DurEpisodes[,5] <- rowSums(DurEpisodes)
  DurEpisodes[numstates+1,] <- colSums(DurEpisodes)
  dimnames(DurEpisodes) <- dimnames(NumEpisodes)


  return (list(n= nsample,
               ne = sum(ns),
               nt = sum(ns-1),
               types=NumEpisodes,
  			   sojourn=DurEpisodes))
 }
