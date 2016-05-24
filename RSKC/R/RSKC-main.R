RSKC <-
  function(d,ncl,alpha,L1=12,nstart=200,silent=TRUE,scaling=FALSE,correlation = FALSE){
    if (alpha > 1 | alpha < 0) stop("alpha must be between 0 and 1")
    if (!is.null(L1)) if (L1<1) stop("L1 value must be greater or equal to 1 or NULL!")
    if (is.data.frame(d)) d <- as.matrix(d)
    r.ncl <- round(ncl)
    if (ncl != r.ncl) ncl <- r.ncl

    if (ncl <= 1) stop("ncl must be positive integer > 1! but ncl=",ncl)

    if (scaling) d=scale(d)
    if (correlation) d = t(scale(t(d)))
    if (is.null(L1)) sparse<-FALSE else{ sparse<-TRUE}
    n<-nrow(d);Nout<-floor(n*alpha)
    f<-ncol(d);g<-f+1
    W<-rep(1,f);sumW<-f # for non-sparse 
    if( sum(is.na(d))==0 ) ## is.na(d) == TRUE if d is na or d is nan
      {
        
        miss<-FALSE
        if(sparse){
          
          Result<-RSKC.a1.a2.b(d,L1,ncl,nstart,alpha,n,f,g,Nout,silent)

        }else{
          ## non sparse
          Result<-RSKC.trimkmeans(d,ncl,trim=alpha,runs=nstart,maxit=10000)
          ## if(Nout!=0) {Result$labels<-class.trimk(d,mu=Result$means,trimC=Result$classification,
          ##     ncl=ncl,Nout=Nout)}
        }
        
      }else{
        d[is.nan(d) ] <- NA
        miss<-TRUE
        if (sparse){
          Result<-RSKC.a1.a2.b.missing(d,L1,ncl,nstart,alpha,n,f,g,Nout,silent)
        }else{
          ## non sparse
          Result<-RSKC.trimkmeans.missing(d=d,ncl=ncl,w=W,trim=alpha,runs=nstart,points=Inf,maxit=10000)
          
          ##if(Nout!=0) {
          ##  Result$labels<-class.trimk.missing(d,mu=Result$means,trimC=Result$classification,
          ##                                     ncl=ncl,Nout=Nout,w=W,sumW=sumW)}
        }
      }         
    ## reported results  
    if(sparse)
      {
        ##sparse
        Result$oW<-sort(Result$oW)
        if(Nout==0){
                                        # sparse K-means
          Result$oW<-Result$oE<-"undefined"
        }
      }else{
                                        #nonsparse
        Result<-modified.result.nonsparse(Result,ncl,f)
        if(Nout==0){
                                        #kmeans
          Result<-modified.result.kmean(Result)
        }
      }   
    Result$disttom<-Result$ropt<-Result$trim<-Result$scaling<-Result$centers<-
      Result$criterion<-Result$classification<-Result$means<-Result$ropt<-Result$k<-Result$runs<-NULL
    
    if (!is.null(colnames(d))) names(Result$weights) <- colnames(d)

    
    Input<-list(N=n,p=f,ncl=ncl,L1=L1,nstart=nstart,alpha=alpha,
                scaling=scaling,correlation=correlation,missing=miss)          
    r2<-c(Input,Result)
    class(r2)<-"rskc"
    return(r2)
  }


modified.result.nonsparse<-function(Result,ncl,f){
  Result$centers<-Result$means;
  Result$oW<-which(Result$classification==ncl+1)
  ##Result$WWSS<-Result$criterion;
  Result$oE<-"undefined";
  Result$weights<-rep(1,f)
  return(Result)
}
modified.result.kmean<-function(Result){
  Result$oE<-Result$oW<-"undefined"
  Result$labels<-Result$classification
  return(Result)
}


## temp <- function(){
##   .Call("tempC",package="RSKC")
## }




