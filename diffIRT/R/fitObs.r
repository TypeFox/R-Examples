fitObs=function(data,margins){
  # this function calculates pr from Eq 3 Maydeu - Olivares & Joe, 2005 (comparable to pi_dot from the matrix at the upper right on p 1010)
  # missing data "NA" allowed
  nit=ncol(data)
  c=combn(nit,margins)
  tmp=(matrix(,1,ncol(c)))
  nms=rep(0,ncol(c))
  base=10^(margins)
  if(margins==1){                                   # gives first moments (proportion correct)
    tmp[1,]=apply(data,2,sum,na.rm=TRUE)
  }
  if(margins>1){                                    # gives higher moments
  for(ii in 1:ncol(c)){
   tmp[,ii]=length(which(apply(data[,c[,ii]],1,sum,na.rm=TRUE)==margins))
  }
  }
return(tmp/nrow(data))
}
