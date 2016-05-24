#' Get order of sequence based on P-value
#' 
#' Takes a sequence as input and find P-value for diffrent orders
#' @param seq - A sequence whose order to be determined
#' @return Returns nothing but prints order of given sequence according to P-value
#' @examples ## Check a first order sequence
#' seq <- c(1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2,1,2) 
#' getOrderPval(seq)
#' 
#' ## Check for second order sequence
#' seq <- c(1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2) 
#' getOrderPval(seq)
#' 
#' ## Check for random order sequence
#' seq <- sample(1:2,50,replace=TRUE)
#' getOrderPval(seq)
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
getOrderPval<- function(seq) {
  n<-length(seq)
  n.states<-length(unique(seq))
  
  ## Preprocessing 
  
  # Determine transition frequency matrix ij
  zero.n<-table(seq)
  first.n<-table(seq[-n],seq[-1]) 
  second.n<-table(seq[1:(n-2)], seq[2:(n-1)], seq[3:n])
  third.n<-table(seq[1:(n-3)], seq[2:(n-2)], seq[3:(n-1)], seq[4:n])
  fourth.n<-table(seq[1:(n-4)], seq[2:(n-3)], seq[3:(n-2)], seq[4:(n-1)], seq[5:n]) 
  fifth.n<-table(seq[1:(n-5)], seq[2:(n-4)], seq[3:(n-3)], seq[4:(n-2)], seq[5:(n-1)], seq[6:n])
  sixth.n<-table(seq[1:(n-6)], seq[2:(n-5)], seq[3:(n-4)], seq[4:(n-3)], seq[5:(n-2)], seq[6:(n-1)], seq[7:n])
  
  # Determine transition frequency matrix for jk 
  first.n.jk <-table(seq[2:(n-1)], seq[3:n])
  second.n.jk <-table(seq[2:(n-2)], seq[3:(n-1)], seq[4:n])
  third.n.jk <-table(seq[2:(n-3)], seq[3:(n-2)], seq[4:(n-1)], seq[5:(n)])
  fourth.n.jk <-table(seq[2:(n-4)], seq[3:(n-3)], seq[4:(n-2)], seq[5:(n-1)], seq[6:(n)])
  fifth.n.jk <-table(seq[2:(n-5)], seq[3:(n-4)], seq[4:(n-3)], seq[5:(n-2)], seq[6:(n-1)], seq[7:(n)])
  
  ##############################################
  # Get p-value for 1st order
  
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      term = (first.n[i,j]) * log(first.n[i,j]/ (zero.n[i] * (zero.n[j] /n ) ) )
      Gsqr = Gsqr + ifelse(is.nan(term),0,term)
    }
  }
  
  Gsqr0 <- (2*Gsqr)
  highOrder <- 1
  lowOrder <- 0
  # Calculate degrees of freedom
  df0 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue01 <- pchisq(Gsqr0,df0,lower.tail=FALSE)
  
  ##############################################
  # Get p-value for 2nd order
  
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        term = (second.n[i,j,k]) * log(second.n[i,j,k]/ (first.n[i,j] * (first.n.jk[j,k] /zero.n[j] ) ) )
        Gsqr = Gsqr + ifelse(is.nan(term),0,term)
      }
    }
  }
  
  Gsqr1 <- (2*Gsqr)
  highOrder <- 2
  lowOrder <- 1
  # Calculate degrees of freedom
  df1 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue12 <- pchisq(Gsqr1,df1,lower.tail=FALSE)
  
  #########################################
  # Get p-value for 3rd order
  
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          term = (third.n[i,j,k,l]) * log(third.n[i,j,k,l]/ (second.n[i,j,k] * (second.n.jk[j,k,l] /first.n[j,k] ) ) )
          Gsqr = Gsqr + ifelse(is.nan(term),0,term)
        }
      }
    }
  }
  Gsqr2 <- (2*Gsqr)
  highOrder <- 3
  lowOrder <- 2
  # Calculate degrees of freedom
  df2 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue23 <- pchisq(Gsqr2,df2,lower.tail=FALSE)
  
  ############################################
  # Get p-value for 4th order
  
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
  Gsqr3 <- (2*Gsqr)
  highOrder <- 4
  lowOrder <- 3
  # Calculate degrees of freedom
  df3 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue34 <- pchisq(Gsqr3,df3,lower.tail=FALSE)
  
  #############################################
  # Get p-value for 5th order
  
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            for ( x in 1:n.states){
              term = (fifth.n[i,j,k,l,m,x]) * log(fifth.n[i,j,k,l,m,x]/ (fourth.n[i,j,k,l,m] * (fourth.n.jk[j,k,l,m,x] /third.n[j,k,l,m] ) ) )
              Gsqr = Gsqr + ifelse(is.nan(term),0,term)
            }
          }
        }
      }
    }
  }
  Gsqr4 <- (2*Gsqr)
  highOrder <- 5
  lowOrder <- 4
  # Calculate degrees of freedom
  df4 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue45<- pchisq(Gsqr4,df4,lower.tail=FALSE)
  
  #########################################
  # Get p-value for 6th order
  
  Gsqr <- 0
  for ( i in 1:n.states){
    for ( j in 1:n.states){
      for ( k in 1:n.states){
        for ( l in 1:n.states){
          for ( m in 1:n.states){
            for ( x in 1:n.states){
              for ( y in 1:n.states){
                term = (sixth.n[i,j,k,l,m,x,y]) * log(sixth.n[i,j,k,l,m,x,y]/ (fifth.n[i,j,k,l,m,x] * (fifth.n.jk[j,k,l,m,x,y] /fourth.n[j,k,l,m,x] ) ) )
                Gsqr = Gsqr + ifelse(is.nan(term),0,term)
              }
            }
          }
        }
      }
    }
  }
  Gsqr5 <- (2*Gsqr)
  highOrder <- 6
  lowOrder <- 5
  # Calculate degrees of freedom
  df5 <- (n.states^highOrder - n.states^lowOrder)*(n.states - 1)
  # Calculate p-value
  Pvalue56<- pchisq(Gsqr5,df5,lower.tail=FALSE)
  
  
  #############################################
  
  results <- matrix(nrow=2,ncol=6,NA)
  results[1,] <- c("0","1","2","3","4","5")
  rownames(results) <- c("Order","P-values")
  colnames(results) <- c("test1","test2","test3","test4","test5","test6")
  results[2,] <- c(Pvalue01,Pvalue12,Pvalue23,Pvalue34,Pvalue45,Pvalue56)
  print(results)
  
  ord = 0
  
  if(Pvalue56 < 0.05){
    print("The given sequence is sixth or higher order")
  } else if (Pvalue45 < 0.05){
    print("The given sequence is of Fifth Order")
  } else if (Pvalue34 < 0.05){
    print("The given sequence is of Fourth Order")
  } else if (Pvalue23 < 0.05){
    print("The given sequence is of Third Order")
  } else if (Pvalue12 < 0.05){
    print("The given sequence is of Second Order")
  } else if (Pvalue01 < 0.05){
    print("The given sequence is of First order")
  } else{
    print('The given sequence is of Random Order')
  }
  
}