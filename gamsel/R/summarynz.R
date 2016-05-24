summarynz=function(object, EPS=0){
  degrees=object$degrees
  p=length(degrees)
  counter=rep(1:p,degrees)
  counter=diag(p)[counter,]
  numlin=apply(abs(object$alpha)>EPS,2, sum)
  anyb=t(counter)%*%abs(object$beta)
  numnlin=apply(anyb>EPS,2,sum)
  numnz=apply( (abs(object$alpha)>EPS)|(anyb>EPS),2,sum)
  cbind(Linear=numlin,Nonlinear=numnlin,Nonzero=numnz)
}

