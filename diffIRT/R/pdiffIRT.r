pdiffIRT=function(y,ai,vi,ter,sd2A,sd2V,A,W,model,eps){
ou=c()
  tmp=function(x) mdiffIRT(x,ai,vi,ter,sd2A,sd2V,A,W,model,eps)
  for(ii in 1:length(y)) 
  if(ii==1) ou[ii]=integrate(tmp,exp(ter)+.001,y[ii])$value
  else ou[ii]=ou[ii-1]+integrate(tmp,y[ii-1],y[ii])$value
return(ou)
}

