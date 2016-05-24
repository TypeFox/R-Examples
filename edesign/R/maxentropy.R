maxentropy <- function(A,nf,ns,method="d",S.start=NULL,rtol=1e-6,
                       mattest=TRUE,etol=0,verbose=FALSE)
  {
    A <- as.matrix(A)
    na <- dim(A)[1]
    nf <- as.integer(nf)
    ns <- as.integer(ns)
    ne <- na-nf
    verbose<-verbose*1

    if (nf==0) stop("algorithm fails for nf=0")
    if (ns==0) stop("number of points have to be added is 0")
    if (ne==0) stop("number of eligible points is 0")

    if(mattest)
      {
        if(!is.double(eigen(A)$values))
          stop("Matrix isn't symmetric")
        if(sum(eigen(A)$values>etol)<na)
          stop("Matrix isn't positive definite")
      }
    
    testmethod <- switch(method, d="ok",g="ok",dc="ok",gc="ok",c="ok","error")
    if(testmethod=="error")
      stop("wrong argument for \"method\", should be one of \"d\", \"g\", \"dc\", \"gc\", \"c\"!")
    
     if(!is.null(S.start))
       {
         if(length(S.start)!=ne)
           {
             stop("Error in dimension of S.start")
           }
         if(length(S.start[S.start!=0])!=ns)
           {
             stop("S.start is infeasible")
           }
       }

    
    if(is.null(S.start))
      {
        ans <- switch(method,
                      d=dualgreedy(A,nf,ns,mattest=FALSE),
                      g=greedy(A,nf,ns,mattest=FALSE),
                      gc=interchange(A,nf,ns,dualgreedy(A,nf,ns,mattest=FALSE)$S),
                      dc=interchange(A,nf,ns,greedy(A,nf,ns,mattest=FALSE)$S)
                      )               
        
        S.start <- ans$S
        det.start <- ans$det
      }
    else
      {
        if(method=="c")
          {
            ans <- interchange(A,nf,ns,S.start)
            S.start <- ans$S
            det.start <- ans$det
          }
        else
          {
            if(nf>0)
              ind <- c(1:nf,S.start[S.start!=0])
            else
              ind <- c(S.start[S.start!=0])
            ni<-length(ind)
            ans.det <- .Fortran("subde1",
                                det=double(1),
                                as.double(A),   #A
                                as.integer(na), #LDA
                                as.integer(na), #NA
                                double(na*na),  #AS
                                as.integer(na), #LDAS
                                as.integer(ind),#IND
                                as.integer(ni), #NI
                                double(na*na),  #BS
                                as.integer(na), #LDBS
                                integer(1),     #IERR
                                PACKAGE="edesign"
                                )
            det.start <- ans.det$det
          }
      }

    rtol <- rtol * det.start

    ans.entrp <- .C("entrp",
                    A=as.double(A),
                    lda=as.integer(na),
                    na=as.integer(na),
                    nf=as.integer(nf),
                    ne=as.integer(ne),
                    ns=as.integer(ns),
                    S=as.integer(S.start),
                    opt=as.double(det.start),
                    integer(ne), # S_Work
                    integer(ne), # F
                    integer(ne), # E
                    integer(na), # ind
                    integer(na), # ind1
                    double(na*na), # As
                    ldas=as.integer(na),
                    double(na*na), # Bs
                    ldbs=as.integer(na),
                    double(na*na), # Cs
                    ldcs=as.integer(na),
                    double(na*na), # Inv
                    ldinv=as.integer(na),
                    double(ns), # W
                    double(8*na), # WORK
                    as.integer(8*na), # LWORK
                    integer(5*na), # IWORK
                    tol=as.double(rtol),
                    maxcount=integer(1),
                    iter=integer(1),
                    verbose=as.integer(verbose),
                    PACKAGE="edesign"
                    )         
    ret <- list(S.start=S.start,det.start=det.start,S=ans.entrp$S,
         det=ans.entrp$opt,maxcount=ans.entrp$maxcount,iter=ans.entrp$iter)
    ret$A <- A
    ret$na <- na
    ret$nf <- nf
    ret$ns <- ns
    ret$ne <- ne
        ret$method <- "maxentropy (exact branch and bound solution)"
    class(ret) <- "monet"
    ret
  }

