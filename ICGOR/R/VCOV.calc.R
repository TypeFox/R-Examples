VCOV.calc <-
function(object){
   P<-length(object$Beta)
   L<-length(object$gl)
   N<-nrow(object$sdata)

hess<-99
if(prod(object$gl>0)==1){
  try1<-try(solve(-object$Hessian),silent=TRUE)
  if(!is.character(try1) & prod(is.finite(try1))==1){
  vcov.bg<-try1[1:P,1:P]
  hess<-0
  }else{
     try2<-try(qr.solve(-object$Hessian),silent=TRUE)
     if(!is.character(try2) & prod(is.finite(try2))==1){
       vcov.bg<-try2[1:P,1:P]
       hess<-1
     }
  }
}
if(hess==99){
  		v<--object$Hessian
  		A = v[1:P, 1:P]
  		B = v[1:P, (P + 1):(P + L)]
  		C = v[(P + 1):(P + L), 1:P]
  		D = v[(P + 1):(P + L), (P + 1):(P + L)]
  		vcov.bg = try(ginv(A - B %*% ginv(D) %*% C),silent=TRUE)
            if(!is.character(vcov.bg) & prod(is.finite(vcov.bg))==1) hess<-2
          }
if(hess==99){
	x0<-c(object$Beta,object$gl)
	Hessian<-try(hessian(f=loglik,x0=x0,sdata=object$sdata,Xp=object$Xp,r=object$r,bl.Li=object$bl.Li,bl.Ri=object$bl.Ri),silent=TRUE)
	VCOV<-diag(solve(-Hessian))
	numer.bg<-try(solve(-Hessian),silent=TRUE)
      if(!is.character(numer.bg)){vcov.bg=diag(numer.bg)[1:P]; hess<-3}
	}
return(list(vcov.bg=vcov.bg,hess=hess))
}
