setClass(
  Class="CURobj", 
  representation=representation(
    C='matrix',
    U='matrix',
    R='matrix',
    C.leverage.score='numeric',
    R.leverage.score='numeric',
    C.index='integer',
    R.index='integer',
    Error='numeric'
  )
)





CUR <- function(A, c=dim(A)[2], r=dim(A)[1], k=NULL, sv=NULL,
                         method="random", alpha=1, weighted=FALSE, beta=4,
                         matrix.return=TRUE, error.return=FALSE)
{

  #function to compute the leverage scores
  levscores <- function(v,k) {
    if (k==1) {
      v[,1]^2
    } else {
      if (weighted) {qq<-sv$d[1:k]^beta; scores=apply(apply(v[,1:k]^2,1,function(x) x*qq),2,sum); scores/sum(scores)} else apply(v[,1:k]^2,1,sum)/k
    }
  }


  method=match.arg(method,c("random","exact.num.random","top.scores","ortho.top.scores","highest.ranks"))
  m=dim(A)[1];n=dim(A)[2]
  ScoresC=NULL
  ScoresR=NULL
  Error=NULL

  if(is.null(k)) {
    if(is.null(sv)) sv=svd(A)
    cs=cumsum(sv$d)
    k=length(cs[cs<cs[length(cs)]*0.8])
  }
  if(is.null(sv)) sv=svd(A,k,k)

  if(is.null(rownames(A))){
    rownames(A)=1:dim(A)[1]
  }
  if(is.null(colnames(A))){
    colnames(A)=1:dim(A)[2]
  }

  if (method=="random") {
    #select columns randomly
    if (c==n) IndC = 1:c else {ScoresC=levscores(sv$v,k); repeat { IndC = which(c*ScoresC>=runif(n)); if (length(IndC)>0) break }} #we repeat until at least 1 column is selected (it is a concern if c is small)
    if (r==m) IndR = 1:r else {ScoresR=levscores(sv$u,k); repeat { IndR = which(r*ScoresR>=runif(m)); if (length(IndR)>0) break }}

  } else if (method=="exact.num.random") {
    #select randomly, but the threshold is adjusted, to pick exatly c columns
    if (c==n) IndC = 1:c else {ScoresC=levscores(sv$v,k); IndC = order(c*ScoresC-runif(n),decreasing=TRUE)[1:c]}
    if (r==m) IndR = 1:r else {ScoresR=levscores(sv$u,k); IndR = order(r*ScoresR-runif(m),decreasing=TRUE)[1:r]}

  } else if (method=="top.scores") {
    #select the columns with highest leverage score
    if (c==n) IndC = 1:c else {ScoresC=levscores(sv$v,k); IndC = order(ScoresC,decreasing=TRUE)[1:c]}
    if (r==m) IndR = 1:r else {ScoresR=levscores(sv$u,k); IndR = order(ScoresR,decreasing=TRUE)[1:r]}

  } else if (method=="ortho.top.scores") {
    if (c==n) {
      IndC = 1:c
    } else {
      pi=levscores(sv$v,k)
      ScoresC=pi
      IndC=integer();
      vnormsqr=apply(A,2,function(x) sum(x^2)) #Norm square of the columns of A: <c_i|c_i>.
      vort=A  #Orthogonal part of the vectors. Initially the whole vector is orthogonal.
      vortnormsqr=vnormsqr #Norm square of the orthogonal part.
      for (i in 1:c)
      {
        IndC[i]=which.max(pi+alpha*sqrt(vortnormsqr/vnormsqr))   #The next column is selected
        sn = vortnormsqr[IndC[i]]    # sn = <v_S|v_S>
        if (sn>0) {
          delta = (vort[,IndC[i]] %*% vort);    #Scalar product of the selected vector with all of the vectors: <v_i|v_S>
          vort = vort - vort[,IndC[i]] %*% delta /sn    #The part, parallel with the new vector is substracted: v_i -> v_i-<v_i|v_S>*v_S/<v_S|v_S>
          vortnormsqr = pmax(vortnormsqr-delta^2/sn,0)    #The norm square is adjusted accordingly: <v_i|v_i> -> <v_i|v_i>-<v_i|v_S>^2/<v_S|v_S>. Rounding errors may set vortnormsqr negative, pmax mitigates this.
        } #(sn>0)
        pi[IndC[i]] = 0   #To avoid being selected again.
      }
    } # c<n
    #and the rows
    if (r==m) {
      IndR = 1:r
    } else {
      pi=levscores(sv$u,k)
      ScoresR=pi
      IndR=integer();
      vnormsqr=apply(A,1,function(x) sum(x^2))
      vort=A  #orthogonal part of the vectors
      vortnormsqr=vnormsqr
      for (i in 1:r)
      {
        IndR[i]=which.max(pi+alpha*sqrt(vortnormsqr/vnormsqr))
        sn = vortnormsqr[IndR[i]]
        if (sn>0) {
          delta = (vort %*% vort[IndR[i],]);
          vort = vort - delta %*% vort[IndR[i],] /sn
          vortnormsqr = pmax(vortnormsqr-delta^2/sn,0)
        } #(sn>0)
        pi[IndR[i]] = 0 #To avoid being selected again.
      }
    } #r<m

  } else if (method=="highest.ranks") {
    if (k==1) stop("For k=1 use method=top.scores instead.")
    if (c==n) {
      IndC = 1:c
    } else {
      crk=apply(apply(sv$v[,1:k]^2,1,cumsum),1,rank) #cumulative rank of the leverage scores
      w=apply(crk,1,which.max) #position of maxima
      mx=apply(crk,1,max) #value of maxima
      ScoresC=w/(n+1)-mx
      IndC=order(-mx,w)[1:c]
      crk=NULL #to release memory
    }
    #and the rows
    if (r==m) {
      IndR = 1:r
    } else {
      crk=apply(apply(sv$u[,1:k]^2,1,cumsum),1,rank) #cumulative rank of the leverage scores
      w=apply(crk,1,which.max) #position of maxima
      mx=apply(crk,1,max) #value of maxima
      ScoresR=w/(m+1)-mx
      IndR=order(-mx,w)[1:r]
      crk=NULL #to release memory
    }

  } else stop("Unknown method: ",method)

#computation of U
  if(matrix.return) {
    if (length(IndC)==1) {
      C = matrix(A[,IndC],ncol=1) #submatrix must be recoerced from vector
    } else {
      C = A[,IndC];    # Choose c' columns from A. 
    }
    if (length(IndR)==1) {
      R = matrix(A[IndR,],nrow=1) #submatrix must be recoerced from vector
    } else {
      R = A[IndR,]      ;    # Choose r' rows from A.
    }
    # Compute U.
    if (c==n) {
      U = MASS::ginv(R)
    } else {         
      if (r==m) {
        U = MASS::ginv(C)
      } else {
        U = MASS::ginv(C) %*% A %*% MASS::ginv(R) ;
      }
    }
    if (error.return) {
      if (r<c) {
        Error=Matrix::norm(A - C %*% U %*% R,"F")
      } else {
        Error=Matrix::norm(A - C %*% (U %*% R),"F")
      }
    }
  } else {
    U=NULL
    C=NULL
    R=NULL
  }

  res = new(Class='CURobj')
  #result=list(C=C, U=U, R=R, IndC=IndC, IndR=IndR, ScoresC=ScoresC, ScoresR=ScoresR, Error=Error)
  if (!is.null(C)) res@C = C
  if (!is.null(ScoresC)) res@C.leverage.score = ScoresC
  if (!is.null(IndC)) res@C.index =IndC
  if (!is.null(R)) res@R = R
  if (!is.null(ScoresR)) res@R.leverage.score = ScoresR
  if (!is.null(IndR)) res@R.index =IndR
  if (!is.null(U)) res@U = U
  if (!is.null(Error)) res@Error = Error
  return(res)
}


setGeneric(
 name= "leverage",
 def=function(object, C=TRUE){
  standardGeneric("leverage")
 }
)

setMethod(
  f="leverage",
  signature= "CURobj",
  definition=function(object, C=TRUE){
    if (C==FALSE){
      res = object@R.leverage.score
    } else if (C==TRUE){
      res = object@C.leverage.score
    }
    return(res)
  }
)


setGeneric(
  name= "plotLeverage",
  def=function(x, C=TRUE, mplr=1000, top.n=100, top.col='red', top.pch=16, ul.col='red', ul.lty=2, ...){
    standardGeneric("plotLeverage")
  }
)

setMethod(
  f="plotLeverage",
  signature= "CURobj",
  definition=function(x, C=TRUE, mplr=1000, top.n=100, top.col='red', top.pch=16, ul.col='red', ul.lty=2, ...){    
    if (mplr!=1 && mplr!=0){
      lab=paste('Leverage score (*', mplr, ')', sep='')
      l.s = leverage(x, C=C)*mplr
    } else {
      lab = 'Leverage score'
      l.s = leverage(x, C=C)
    }
    plot(l.s, type='h', ylab=lab, ...)
    idx = rev(sort.list(l.s))[1:top.n]
    points(idx, l.s[idx], pch=top.pch, col=top.col)
    abline(h=mean(l.s), col=ul.col, lty=ul.lty)
  }
)

setGeneric(
  name= "topLeverage",
  def=function(object, C=TRUE, top.n=100, sort=TRUE){
    standardGeneric("topLeverage")
  }
)

setMethod(
  f="topLeverage",
  signature= "CURobj",
  definition=function(object, C=TRUE, top.n=100, sort=TRUE){
    l.s = leverage(object, C=C)
    idx = rev(sort.list(l.s))[1:top.n]
    if (sort==TRUE){
      idx = sort(idx)
    }
    return(idx)
  }
)

setGeneric(
  name= "getR",
  def=function(object){
    standardGeneric("getR")
  }
)

setMethod(
  f="getR",
  signature= "CURobj",
  definition=function(object){
    return(object@R)
  }
)

setGeneric(
  name= "getU",
  def=function(object){
    standardGeneric("getU")
  }
)

setMethod(
  f="getU",
  signature= "CURobj",
  definition=function(object){
    return(object@U)
  }
)

setGeneric(
  name= "getC",
  def=function(object){
    standardGeneric("getC")
  }
)

setMethod(
  f="getC",
  signature= "CURobj",
  definition=function(object){
    return(object@C)
  }
)

setGeneric(
  name= "getError",
  def=function(object){
    standardGeneric("getError")
  }
)

setMethod(
  f="getError",
  signature= "CURobj",
  definition=function(object){
    return(object@Error)
  }
)

