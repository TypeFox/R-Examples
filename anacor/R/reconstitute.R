`reconstitute` <-
function(tab,order=0,eps=1e-6,verbose=FALSE)
{
# NORA's algorithm for missings

  x0 <- tab
  x0[is.na(x0)] <- 0                                    #replace NA with 0's		
  xold<-x0                                              #save original table
  fold<- -Inf
  repeat {
        rx<-rowSums(xold)                               #row margins new table
        cx<-colSums(xold)                               #column margins new table
        tx<-sum(xold)                                   #grand total new table
  	y<-outer(rx,cx)/tx                              #frequencies under independence y

        xnew <- tab
        xnew[is.na(tab)] <- y[is.na(tab)]               #replace NA's by y
  	fnew <- sum(x0*log(ifelse(y==0,1,y)))           #

        if (verbose) print(c(fold,fnew))

        if ((fnew - fold) < eps) break()
  	xold<-xnew                                      #update new table
        fold<-fnew
  }
  return(xnew)
}

