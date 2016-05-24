

htvar.list<-function(xcheck, Dcheck){
    rval<-sapply(Dcheck, function(stagei)
             {htvar.matrix(rowsum(xcheck,stagei$id),stagei$dcheck)})
    rval
}

## used in twophase2var()
htvar.matrix<-function(xcheck, Dcheck){
  if (is.null(dim(xcheck)))
    xcheck<-as.matrix(xcheck)
  rval<-apply(xcheck,2, function(xicheck)
              apply(xcheck,2, function(xjcheck)
                    as.matrix(Matrix::crossprod(xicheck, Dcheck%*%xjcheck))
                    ))
  if(is.null(dim(rval))) dim(rval)<-c(1,1)
  rval
}

## used in ppsvar, twophase2var
ygvar.matrix<-function(xcheck,Dcheck){
  ht<-htvar.matrix(xcheck,Dcheck)
  if (is.null(dim(xcheck))){
    corr <- sum(Dcheck%*%(xcheck*xcheck))    
  } else {
    corr <- apply(xcheck,2, function(xicheck)
                apply(xcheck,2, function(xjcheck)
                      sum(Dcheck%*%(xicheck*xjcheck))
                      ))
  }
  rval<-ht-corr 
}


