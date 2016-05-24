"ftf.x.fit" <-function(y,x,S,control){
  n=dim(S)[1]
  m=length(y)

  L=1:m
  U=(m+1):n
    
  V=diag(n)-S
  Q0=V%*%x
  Q1=t(x)%*%V
  cv=try(solve(t(x)%*%V%*%x,Q1),silent=TRUE)
  if(class(cv)=="try-error"){
    if(control$warn)
      warning("Matrix Inverse Problem ... Using ridge with l=0.1\n")
    cv=try(solve(t(x)%*%V%*%x+0.1*diag(dim(x)[2]),Q1),silent=TRUE)
    if(class(cv)=="try-error"){
      if(control$warn)
        warning("Matrix Inverse Problem ... Using ginv instead of solve\n")
      cv=ginv(t(x)%*%V%*%x)%*%Q1
    }
  }
  
  H=Q0%*%cv+S
  if( n==m){
    tdf=sum(diag(H))
    yhat=as.vector(H%*%y)
    coefs=cv%*%y
    fg=yhat-x%*%coefs
    xstr=list(tdf=tdf,yhat=yhat,coefs=coefs,fg=fg) 
   return(xstr)
  }
  if( (n-m)==1){
    HUUi=1/(1-H[U,U])
    A=H[L,L]+HUUi*outer(H[L,U],H[U,L])
    tdf=sum(diag(A))
    yhat=c(A%*%y[L],HUUi*H[U,L]%*%y[L])
    ystar=c(y[L],yhat[U])
    coefs=cv%*%ystar
    fg=yhat-x%*%coefs
    xstr=list(tdf=tdf,yhat=yhat,coefs=coefs,fg=fg)
    return(xstr)
  }
  HUUi=try(solve(diag(n-m)-H[U,U]),silent=TRUE)
  if(class(HUUi)=="try-error"){
    if(control$warn)
      warning("Matrix Inverse Problem ... Using ginv instead of solve\n")
    HUUi=ginv(diag(n-m)-H[U,U])
  }
  B=HUUi%*%H[U,L]
  A=H[L,L]+H[L,U]%*%B
  tdf=sum(diag(A))
  yhat=as.vector(c(A%*%y,B%*%y))
  ystar=c(y[L],yhat[U])
  coefs=cv%*%ystar
  fg=yhat-x%*%coefs
  xstr=list(tdf=tdf,yhat=yhat,coefs=coefs,fg=fg)
  return(xstr)
}
