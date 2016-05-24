#' Determine if Time Homogeneity is present in the give sequence.
#' 
#' Takes a sequence as input and finds if Time Homogeneity is present or not.
#' @param seq - A sequence whose Time Homogeneity is to be determined
#' @return Returns nothing but prints output representing presence or absence of Time Homogeneity
#' @examples ## Check for a homogenous sequence
#' seq <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2) 
#' checkTH(seq)
#' 
#' ## Check for a heterogenous sequence
#' seq <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1)
#' checkTH(seq)
#' @references
#' [1] Markov Chain Test for Time Dependence and Homogeneity: An Analytical and Empirical Evaluation 
#' Baris Tan and Kamil Yilmaz European Journal of Operational Research 137 (2002) 524-543
#' @export
checkTH <- function(seq) {
  
  # Preprocessing
  n<-length(seq)  # Get length
  n.states<-length(unique(seq))   # Get number of states
  
  zero.n<-table(seq)  # Transition frequency matrix, 0th order
  zero.p<-zero.n/n  # Transition probability matrix. 0th order
  
  first.n<-table(seq[-n],seq[-1])
  first.p<-first.n/rowSums(first.n)
  
  # For checking Time Homogeniety we would divide given sequence into multiple segments 
  # and compate properties of each segment with others. 
  
  # At first we will divide the sequence in 3 segments then 4 and then 5 to confirm our 
  # results are independent of number of segments we create from given sequence
  
  segCountV <- c(3,4,5)  
  
  # So first we will divide sequence in 3 segments, for next iterations in 4 and then in 5 segments.
  
  # To store results from each iteration about TH test
  resVector <- c(0,0,0)   # By default we assume TH not present
  
  ### Check Homogeneity iteratively  #####################################################
  for (itr in 1:3){
    # Get number of segments for this iteration
    segCount <- segCountV[itr] 
    # Calculate segment length for this iteration
    segLen <- floor(n/segCount) 
    # In SegMat (Segment Matrix) each column is a segment of original sequence
    segMat <- matrix(nrow=segLen,ncol=segCount,NA) 
    
    # Get values in SegMat
    for (i in 1:segCount){
      low <- (segLen*(i-1))+1
      high <- segLen*i
      segMat[,i] = seq[low:high]
    }
    
    ## We can access i-th segment (column) as segMat[,i][1:segLen]
    
    ## Lets declare transition probability matrix for each segment
    first.p.seg = first.n.seg <- array(0, dim = c(segCount,n.states,n.states))
    
    ## Lets fill in the 'n' matrix for each segment
    for (k in 1:segCount){
      first.n.seg[k,,] <- table(segMat[,k][1:segLen-1],segMat[,k][2:segLen])
    }
    
    ## Lets fill in 'p' for each segment
    for (k in 1:segCount){
      first.p.seg[k,,] <-first.n.seg[k,,]/rowSums(first.n.seg[k,,])
    }
    
    Lambda <- 0
    
    for (k in 1:segCount){
      for (i in 1:n.states){
        for (j in 1:n.states){
          term <- first.n.seg[k,i,j] *(log(first.p.seg[k,i,j]) - log(first.p[i,j]))
          Lambda = Lambda + ifelse(is.nan(term),0,term)
        }
      }
    }
    
    chisq <- 2*Lambda
    # Calculate degrees of freedom
    df <- (segCount - 1)*(n.states -1)*(n.states)
    # Calculate p-value
    pValue <- pchisq(chisq,df,lower.tail=FALSE)
    
    # If p-value is greater than or equal to 0.5 then TH present
    if (pValue >= 0.05){
      resVector[itr] = 1
    }
    
  } ### End of iterative process
  
  
  ## If all iterations indicate the presence of time homogeneity only then we will conclude positively
  
  if(resVector[1]*resVector[2]*resVector[3]){  
    print('Time homogeneity present in given sequence')
  }else{
    print('Time homogeneity absent in given sequence')
  }
}
