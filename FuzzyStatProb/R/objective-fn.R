.unos = function(m){
	temp=m;
	d = dim(m)[2];
	for(j in 1:d)
    temp[j,d]=1
  temp
}
.stationary = function(x,d,listaposx,listaposy,component,lowerboundsmatrix,upperboundsmatrix,maximize){

  P<-matrix(0,d,d);
  acum=0;
  indicex=1;
  for(i in 1:length(listaposx)){
    if(listaposy[i]<0){   # we have ended a row
      P[listaposx[i],(-1)*listaposy[i]] = 1 - acum;
      if((1-acum) < lowerboundsmatrix[listaposx[i],(-1)*listaposy[i]] - 1E-5 ||
         (1-acum) > upperboundsmatrix[listaposx[i],(-1)*listaposy[i]] + 1E-5){     
         if(maximize==TRUE){ 
            return(0); 
         }
         else{
            return(1); 
         }
      }   
      acum = 0;
    }
    else{ 
      P[listaposx[i],listaposy[i]] = x[indicex];
      acum = acum + x[indicex];
      indicex = indicex + 1;
    }
  }
  # ------------------------------------------------------------------------
  #  Computation of stationary probabilities with crisp transition matrix P
  # ------------------------------------------------------------------------
  identidad<-diag(d);
  ceros<-matrix(0,d,d);
  mat2=.unos(ceros) %*% solve(.unos(P-identidad));
  vec=mat2[1,];  # Now vec contains the stationary distribution
  # ------------------------------------------------------------------------
  if(maximize==TRUE){
    return((-1)*vec[component]);
  }
  else{ 
    return(vec[component]);
  }
}