"rsunpack" <-
function(x,n,k,nwords,userfunc,...){
     # check for NAs (-2^31 not defined in R) in simulated matrices
     # set value to 0

     nas<-FALSE
     if (k>=32) {
       idx<-(1:length(x))[is.na(x)]   # indexvector for NAs
       nas<-(length(idx)>0)
       x[idx]<-0
     }

     t<-vector(length=n*k)

     # calls unpacking routine
     out<-.Fortran("unpack",
                 as.integer(x),
                 as.integer(nwords),
                 mat=as.integer(t),
                 as.integer(n),
                 as.integer(k)
     )
     m<-matrix(out$mat,nrow=n)
     # replace NAs with bitpattern corresponding to -2^31,
     # i.e., 0 0 0.... 0 1
     if (nas) {
        idx1 <- ceiling(idx/nwords)                  # index for rows
        targetbyte <- idx%%nwords                    # which byte in row is affected
        last <- k%%32                                # last column in targetbyte
        idx2 <- (targetbyte - 1*(last!=0))*32 + last # index for column
        m[idx1,idx2]<-1
     }
     # calls user function to calculate statistic(s)
     RET<-userfunc(m,...)
     RET
}

