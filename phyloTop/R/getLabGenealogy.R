#' Create genealogy
#' 
#' Create a labelled genealogy from an epidemiological record
#' 
#' @author Caroline Colijn \email{c.colijn@imperial.ac.uk}
#' @author Michelle Kendall \email{michelle.louise.kendall@@gmail.com}
#'   
#' @param epirecord an epidemiological record, as output from the function \code{\link{makeEpiRecord}}. It must be a matrix of at least two rows and with five columns named "Infectee", "Infector", "InfnTime", "RecTime", "DoneFlag".
#' @param epsilon an optional small number to be used for branch lengths which would otherwise be zero.
#'
#' @return An object of class \code{phylo} representing the transmission tree from infectors to infectees.
#' 
#' @seealso \code{\link{makeEpiRecord}}
#' 
#' @examples
#' ## Generate an epidemiological record:
#' myepirecord <- makeEpiRecord(c(1,2,3,4))
#' ## make the corresponding genealogy from this record:
#' mygenealogy <- getLabGenealogy(myepirecord)
#' plot(mygenealogy)
#' 
#' @export
getLabGenealogy <- function(epirecord,epsilon=0.001234) {
  # make genealogy from a matrix with columns "Infectee", "Infector", "InfnTime", "RecTime", "DoneFlag" as created by makeTransTree
  cols <- c("Infectee","Infector","InfnTime","RecTime","DoneFlag")
    
  # check we have expected column names
  if (!setequal(colnames(epirecord),cols)) {
    stop(paste('The object epirecord must have columns named "Infectee", "Infector", "InfnTime", "RecTime", "DoneFlag"', sep=''))
  }
  
  # check epirecord has more than one row
  if (length(epirecord[,1])==1) {
    stop(paste0("The object epirecord has only one row; cannot build a meaningful tree from this information"))
  }
    
  # sort the rows of the matrix epirecord so that InfnTime goes from most recent to oldest (used to call the function sortmyepi but this is quicker)
  sortedepirecord <- epirecord[order(epirecord[,"InfnTime"],decreasing=TRUE),]

  # extract columns of sorted epirecord
  Infectors <- sortedepirecord[,"Infector"];
  Infectees <- sortedepirecord[,"Infectee"];
  InfTimes <- sortedepirecord[,"InfnTime"];
  SortedRecTimes <- sortedepirecord[,"RecTime"]; 
  SortedFLAGS <- sortedepirecord[,"DoneFlag"];
    
  # In the event that any sorted recovery times equal some infection times, add epsilon to these times so that they no longer clash.
  SortedRecTimes[which(SortedRecTimes==InfTimes)] <- SortedRecTimes[which(SortedRecTimes==InfTimes)] + epsilon
  
  NN <- length(InfTimes)-1; # NNodes (internal)
  NumLeaves <- NN+1; # number of tips
  BranchNums <- (NumLeaves+1):(2*NumLeaves-1); # IDs for internal branches of the genealogy
  LengthstoInternals <- rep(0,NN)
  TipLengths <- rep(0,NumLeaves)
  B <- matrix(0,nrow=NN, ncol=2) # use a 2-col [desc desc] format; convert afterwards
    
    
  for (n in 1:NN ) {
    #### PART 1 : connecting the internal branches / tips
    if (n ==1 ) { 
      B[n,]<-c(Infectors[n], Infectees[n]);
    } 
    else { 
      # 1: did the infectEE go on to infect anyone else? 
      Onward <- is.element(Infectees[n],Infectors[1:(n-1)])
        if (Onward) {
          Desc1 <- BranchNums[max(which(Infectors[1:(n-1)] == Infectees[n]))]
        }
        else {
          Desc1 <- Infectees[n]
        }
        
      # 2: did the infectOR go on to infect anyone else? 
      Onward <- is.element(Infectors[n],Infectors[1:(n-1)])
        if (Onward){
          Desc2 <- BranchNums[max(which(Infectors[1:(n-1)]==Infectors[n]))]
        } 
        else {
          Desc2 <- Infectors[n]
        }
      B[n,] <- c(Desc1,Desc2) # if n is NN, then this branch is the root. 
    } 
      
      
    ### PART 2: lengths of internal branches    
    HasInfBefore <- is.element(Infectors[n],Infectors[(n+1):length(Infectors)])
      if (HasInfBefore){
        IND=n+ min(which(Infectors[(n+1):length(Infectors)]==Infectors[n])) # most recently
        LengthstoInternals[n]=InfTimes[n]-InfTimes[IND]; # corresponds to int pt n, 
      } 
      else {
        IND<- which(Infectees==Infectors[n]) # index when this guy got infected
        LengthstoInternals[n]=InfTimes[n]-InfTimes[IND] 
      }
  } # end for loop 
    
  
  ### now assign lengths of tip branches
  for (n in 1:length(Infectees)) {
    HasInfBefore<-is.element(Infectees[n],Infectors)  
    if (HasInfBefore) {
      IND=min(which(Infectors==Infectees[n])) 
      TipLengths[n] <- SortedRecTimes[n]-InfTimes[IND] 
    } 
    else {
      TipLengths[n]<-SortedRecTimes[n]-InfTimes[n]
    }
  } 
  # TipLengths are in the same order as Infectee. We sort the flags along with epirecord so that flag[n] will correspond to TipLength[n]
  ReorderedTipLengths=rep(0,length(TipLengths))
  ReorderedTipLengths[Infectees]=TipLengths    
    

  
  ## now create the genealogy: change format
    
  # initialise
  Edges=matrix(0,nrow=2*NN, ncol=2)
  Lengths<-0*(1:(2*NN))
  FLAGS <- Lengths
    
  ROSortedFLAGS <- SortedFLAGS
  ROSortedFLAGS[Infectees] <- SortedFLAGS

  for (n in 1:NN) {
    Edges[2*n-1,1] <- n+NN+1
    Edges[2*n,1] <- n+NN+1
    Edges[2*n-1,2] <- B[n,1]
    Edges[2*n,2] <- B[n,2]
  }
    
  Lengths[1] <- 0;
  for (n in 1:(2*NN)) { 
    IsDescTip <- (Edges[n,2]<= NN+1)
    if (IsDescTip) {
      Lengths[n]<- ReorderedTipLengths[Edges[n,2]] # if he is a tip, here is his length
      FLAGS[n] <- ROSortedFLAGS[Edges[n,2]] # and here is his flag.  # CONFIRMED
    } 
    else {
      Lengths[n]<- LengthstoInternals[Edges[n,2]-NN-1] 
    }
  } 
  
  # if any lengths are zero, make them length epsilon
  Lengths[Lengths==0]=epsilon
    
  Edges[Edges>NumLeaves]= 3*NumLeaves - Edges[Edges>NumLeaves]; # so branches start at ROOT
  Genealogy <- makePhyloTree(Edges,Lengths,NumLeaves+1) # NOW root is first branch, not last
  return(Genealogy)
}