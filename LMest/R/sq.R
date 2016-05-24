sq <- function(J,t=NULL){

  M=NULL
  if(J == 0) return
  if(!is.null(t)){           # t is present
    if(t == J) M = matrix(1,1,J) else{         # all elements equal to 1
      if(t>1){       # se sono sia 1 che 0
        for(i in 1:(J-t+1)){ 
          S = sq(J-i,t-1)
          r = nrow(S)
          T1 = cbind(matrix(0,r,i-1),rep(1,r),S)
          M = rbind(T1,M)
        }
      }else{
        if(t==1){      # only one element       
          M = matrix(0,J,J)
          for(j in 1:J) M[j,J-j+1] = 1
        }else{ 
          M = matrix(0,1,J) # no elements
        }
      }
    }
  }else{                # there is no t
    if(J==1) M = matrix(c(0,1),2,1) else{
      T1 = sq(J-1);  nt = nrow(T1)
      M = rbind(cbind(rep(0,nt),T1),cbind(rep(1,nt),T1))
    }
  }
  M
}