## methods
seMean<-function(x,...) UseMethod("seMean")
seVar<-function(x,...) UseMethod("seVar")
seDMean<-function(x,...) UseMethod("seDMean")
seDMeanG<-function(x,...) UseMethod("seDMeanG")
seDVar<-function(x,...) UseMethod("seDVar")
seRMean<-function(x,...) UseMethod("seRMean")
seRVar<-function(x,...) UseMethod("seRVar")

 
## do not forget the ... if this is declared in the method!
seMean.default<-function(x,...) {
  sqrt(var(x)/length(x))
}

seVar.default<-function(x,...) {
  sqrt(var((x-mean(x))^2)/length(x))
}

seDMean.default<-function(x,y,rho=1,...) {
  sqrt(var(x)/length(x)+rho^2*var(y)/length(y))
}

seDMeanG.default<-function(x,y,...) {
  nx<-length(x)
  ny<-length(y)
  sqrt(((nx-1)*var(x)+(ny-1)*var(y))/(nx+ny-2)*(1/nx+1/ny))
}
  
seDVar.default<-function(x,y,rho=1,...) {
  sqrt(var((x-mean(x))^2)/length(x)+rho^2*var((y-mean(y))^2)/length(y))
}
  
seRMean.default<-function(x,y,r0,...) {
  if (missing(r0)) r0 <- mean(x)/mean(y)
  sqrt(var(x)/length(x)+r0^2*var(y)/length(y))/abs(mean(y))
}
  
seRVar.default<-function(x,y,r0,...) {
  if (missing(r0)) r0 <- var(x)/var(y)
  sqrt(var((x-mean(x))^2)/length(x)+r0^2*var((y-mean(y))^2)/length(y))/var(y)
}
