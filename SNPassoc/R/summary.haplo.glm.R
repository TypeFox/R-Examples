`summary.haplo.glm` <-
function(object, ...)
 {
   o <- object
   coe<-o$coe
   se<-sqrt(diag(o$var.mat)[1:length(coe)])
   z<-coe/se
   p<-2-2*pnorm(abs(z))
   list(coeficients=cbind(coe,se,z,p))
 }

