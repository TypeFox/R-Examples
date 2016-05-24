VCOV.calc <-
function(object){
   P<-length(object$Beta)
   M<-length(object$Eta)
   L<-length(object$gl)
   N<-nrow(object$sdata)

hess<-99
if(prod(object$gl>0)==1){
  try1<-try(solve(-object$Hessian),silent=TRUE)
  if(!is.character(try1) & prod(is.finite(try1))==1){
  vcov.bg<-try1[1:(M+P),1:(M+P)]
  hess<-0
  }else{
     try2<-try(qr.solve(-object$Hessian),silent=TRUE)
     if(!is.character(try2) & prod(is.finite(try2))==1){
       vcov.bg<-try2[1:(M+P),1:(M+P)]
       hess<-1
     }
   }
}
if(hess==99){
           v<--object$Hessian
           A = v[1:(M+P), 1:(M+P)]
           B = v[1:(M+P), (M+P + 1):(M+P + L)]
           C = v[(M+P + 1):(M+P + L), 1:(M+P)]
           D = v[(M+P + 1):(M+P + L), (M+P + 1):(M+P + L)]
           vcov.bg = try(ginv(A - B %*% ginv(D) %*% C),silent=TRUE)
           if(!is.character(vcov.bg) & prod(is.finite(vcov.bg))==1) hess<-2
          }
if(hess==99){
  x0<-c(object$Eta,object$Beta,object$gl)
  Hessian<-try(hessian(f=loglik,x0=x0,r=object$r,sdata=object$sdata,Xp=object$Xp,Zp=object$Zp,bl.Li=object$bl.Li,bl.Ri=object$bl.Ri),silent=TRUE)
  numer.bg<-try(solve(Hessian)[1:(M+P),1:(M+P)],silent=TRUE)
  if(!is.character(numer.bg)){vcov.bg=numer.bg; hess<-3}
}
return(list(vcov.bg=vcov.bg,hess=hess))
}
