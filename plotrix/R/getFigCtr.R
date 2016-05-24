getFigCtr<-function(pos=c(0.5,0.5)) {
 pars<-par(c("usr","mar","fin","pin"))
 pxspan<-diff(pars$usr[1:2])
 fxspan<-pxspan*pars$fin[1]/pars$pin[1]
 figxctr<-
  pars$usr[1]-(fxspan-pxspan)*pars$mar[2]/(pars$mar[2]+pars$mar[4]) + fxspan*pos[1]
 pyspan<-diff(pars$usr[3:4])
 fyspan<-pyspan*pars$fin[2]/pars$pin[2]
 figyctr<-
  pars$usr[1]-(fyspan-pyspan)*pars$mar[1]/(pars$mar[1]+pars$mar[3]) + fyspan*pos[2]
 return(c(figxctr,figyctr))
}
