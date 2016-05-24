#' Get order of sequence based on BIC value
#' 
#' Takes a sequence as input and find BIC value for diffrent orders
#' @param seq - A sequence whose order to be determined
#' @return Returns nothing but prints order of given sequence according to BIC value
#' @examples ## Check a first order sequence
#' seq <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2) 
#' getOrderBIC(seq)
#' 
#' ## Check for second order sequence
#' seq <- c(1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2) 
#' getOrderBIC(seq)
#' 
#' ## Check for random order sequence
#' seq <- sample(1:2,50,replace=TRUE)
#' getOrderBIC(seq)
#' @references
#' [1] Estimating the order of Markov chain Richard Katz Technometrics, 
#'    vol 12 no 3 (August 1981) pp 243-249 
#' 
#' [2] Determination of the Order of a Markov Chain L.C.Zhao, C.C.Y.Dorea and 
#'    C.R.Goncalves Statistical inference for stochastic processes4, 2001 pp 273-282 
#' 
#' [3] Statistical inference about Markov Chain T.W.Anderson and Leo.A.Goodman. 
#'    The Annals of Mathematical Statistics, Vol 28, No 1 (March 1957), pp89-110
#' @export
getOrderBIC<- function(seq) {
  
  n<-length(seq)
  n.states<-length(unique(seq))
  
  ## Preprocessing 
  
  # Determine transition frequency matrix ij
  zero.n<-table(seq)
  first.n<-table(seq[-n],seq[-1])
  second.n<-table(seq[1:(n-2)], seq[2:(n-1)], seq[3:n])
  third.n<-table(seq[1:(n-3)], seq[2:(n-2)], seq[3:(n-1)], seq[4:n])
  fourth.n<-table(seq[1:(n-4)], seq[2:(n-3)], seq[3:(n-2)], seq[4:(n-1)], seq[5:n])
  
  # Determine transition frequency matrix for jk 
  first.n.jk <-table(seq[2:(n-1)], seq[3:n])
  second.n.jk <-table(seq[2:(n-2)], seq[3:(n-1)], seq[4:n])
  third.n.jk <-table(seq[2:(n-3)], seq[3:(n-2)], seq[4:(n-1)], seq[5:(n)])
  
  # As we consider higher order complexity increase significantly, so we will test up to 3rd order.
  
  # Order check for 0 vs 4
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            term = (fourth.n[i,j,k,l,m]) * log(fourth.n[i,j,k,l,m]/ (third.n[i,j,k,l] * (zero.n[j] /n ) ) )
            Gsqr = Gsqr + ifelse(is.nan(term),0,term)
          }
        }
      }
    }
  }
  
  
  highOrder <- 4
  lowOrder <- 0
  # Calculate degrees of freedom
  df0 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  
  # For 0th order
  BIC0<-(2*Gsqr - (df0)*log(n))
  
  
  # Order check for 1 vs 4
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            term = (fourth.n[i,j,k,l,m]) * log(fourth.n[i,j,k,l,m]/ (third.n[i,j,k,l] * (first.n.jk[j,k] /zero.n[j] ) ) )
            Gsqr = Gsqr + ifelse(is.nan(term),0,term)
          }
        }
      }
    }
  }
  
  
  highOrder <- 4
  lowOrder <- 1
  # Calculate degrees of freedom
  df1 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  
  # For 1st order
  BIC1<-(2*Gsqr - (df1)*log(n))
  
  # Order check for 2 vs 4
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            term = (fourth.n[i,j,k,l,m]) * log(fourth.n[i,j,k,l,m]/ (third.n[i,j,k,l] * (second.n.jk[j,k,l] /first.n[j,k] ) ) )
            Gsqr = Gsqr + ifelse(is.nan(term),0,term)
          }
        }
      }
    }
  }
  
  highOrder <- 4
  lowOrder <- 2
  # Calculate degrees of freedom
  df2 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  
  # For 2nd order
  BIC2<-(2*Gsqr - (df2)*log(n))
  
  # Order check for 3 vs 4
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            term = (fourth.n[i,j,k,l,m]) * log(fourth.n[i,j,k,l,m]/ (third.n[i,j,k,l] * (third.n.jk[j,k,l,m] /second.n[j,k,l] ) ) )
            Gsqr = Gsqr + ifelse(is.nan(term),0,term)
          }
        }
      }
    }
  }
  
  highOrder <- 4
  lowOrder <- 3
  # Calculate degrees of freedom
  df3 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  
  # For 3nd order  
  BIC3<-(2*Gsqr - (df3)*log(n))    
  
  ## Creating a matrix of all AIC values
  
  results <- matrix(nrow=2,ncol=4,NA)
  rownames(results) <- c("Order","BIC")
  colnames(results) <- c("test1","test2","test3","test4")
  results[1,] <- c("0","1","2","3")
  results[2,] <- c(BIC0,BIC1,BIC2,BIC3)
  print(results)
  
  # Printing order based on result matrix
  ord <- which.min(results[2,])-1
  if (ord == 0){
    print("The given sequence is a random sequence")
  } else if(ord == 1 | ord == 2){
    print(paste("The order for the given sequence is: ",ord,sep=''))
  } else {
    print("The sequence is of Third or Higher Order")
  }
  
}