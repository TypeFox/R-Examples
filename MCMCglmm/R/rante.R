"rante"<-function(y,nk=1,n=1,cbeta=TRUE, cvar=TRUE, Ainv=NULL){
 
  if(is.matrix(y)==FALSE){stop("y must be a matrix")}
  if(!nk<ncol(y)){stop("nk must be less than ncol(y)")}
  if(nk<0){stop("nk must be a positive integer")}
  if(abs(nk - round(nk)) > .Machine$double.eps^0.5){stop("nk must be a positive integer")}
  if(n<0){stop("n must be a positive integer")}
  if(abs(n - round(n)) > .Machine$double.eps^0.5){stop("n must be a positive integer")}
  if(!is.logical(cbeta)){stop("cbeta should be TRUE or FALSE")}
  if(!is.logical(cvar)){stop("cbeta should be TRUE or FALSE")}

  if(is.null(Ainv)){
    iAP<-pAP<-xAP<-nzmaxP<-0
  }else{
    pAP<-Ainv@p             
    iAP<-Ainv@i	         
    xAP<-Ainv@x	    
    nzmaxP<-length(xAP)  
  }

  GinvP<-rep(0,n*ncol(y)^2)

  output<-.C("rante",
     as.double(c(y)),  
     as.integer(ncol(y)),          
     as.integer(nrow(y)),          
     as.integer(nk),
     as.integer(n),           
     as.integer(cbeta),
     as.integer(cvar),
     as.integer(iAP),             
     as.integer(pAP),	         
     as.double(xAP),	    
     as.integer(nzmaxP),      
     as.double(GinvP)
  )   
  Ginv<-mcmc(t(matrix(output[[12]],ncol(y)^2,n)))
  return(Ginv)
}
