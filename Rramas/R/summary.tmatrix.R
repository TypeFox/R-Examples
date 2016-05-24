summary.tmatrix <-
function(object,...){
  name.mat<-deparse(substitute(object))
  x <- object
  di <- dim(x)[1]
  m.names <- dimnames(x)[[1]] # ya en as.tmatrix
  # if(is.null(m.names)) m.names <- paste("stage.",1:di ,sep="")
  ea<- eigen(x)
  lambda <-abs( ea$values[1]) # finite rate of increase
  #stable stage/age distribution:
  ssd <- abs(ea$vectors[,1]/sum(ea$vectors[,1]) ) 
  ae <- eigen(t(x))
  vr <- abs(ae$vectors[,1]/ae$vectors[1,1] )# reproductive value
  sensitivity <-  (vr  %*%  t(ssd))  / (t(vr) %*% ssd)[1,1]
  elasticity <- sensitivity * x / lambda
   
  result<- list(lambda=lambda, stable.stage.distribution = ssd,
               reproductive.value =vr, sensitivity = sensitivity,
               elasticity=elasticity,name.mat=name.mat,m.names= m.names)
  class(result)=c("summary.tmatrix", class(result))
  return (result)
 }

