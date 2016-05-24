multinomialCI <- function(x,alpha,verbose=FALSE){
  
  n = sum(x, na.rm=TRUE);
  k = length(x);
  p = x/n;
  c = 0;
  pold=0;
  for(cc in 1:n){
   p = .truncpoi(cc,x,n,k);
   if(p > 1-alpha && pold < 1-alpha) { c = cc; break; };
   pold=p;
  }

  salida = matrix(0,k,2);
  delta=(1-alpha-pold)/(p-pold);
  out=matrix(0,k,5);
  num=matrix(0,k,1);
  c=c-1;  
  vol1=1;
  vol2=1;
  for(i in 1:k){
   num[i,1]=i;
   obsp=x[i]/n;
   out[i,1]=obsp;
   out[i,2]=obsp-c/n;
   out[i,3]=obsp+c/n+2*delta/n;
   if(out[i,2]<0){ out[i,2]=0; }
   if(out[i,3]>1){ out[i,3]=1; }
   out[i,4]=obsp-c/n-1/n;
   out[i,5]=obsp+c/n+1/n;
   if(out[i,2]<0){ out[i,2]=0; }
   if(out[i,3]>1){ out[i,3]=1; }
   vol1=vol1*(out[i,3]-out[i,2]);
   vol2=vol2*(out[i,5]-out[i,4]);
   
   salida[i,1] = out[i,2];
   salida[i,2] = out[i,3];
  }
  c1=c('PROPORTION', 'LOWER(SG)','UPPER(SG)','LOWER(C+1)','UPPER(C+1)');
  cov=100*(1-alpha);
  sg=(x+delta)/n;
  c2=c('SG-midpoint');
  if(verbose==TRUE){
    print('-------------------------------------------------------------');
    print(paste('    ',cov,'% SIMULTANEOUS CONFIDENCE INTERVALS'));
    print('       BASED ON THE METHODS OF SISON AND GLAZ');
    print('-------------------------------------------------------------');
    print(paste('C = ',c));
    print(paste('P(c+1) = ',p));
    print(paste('P(c)   = ',pold));
    print(paste('delta =  ',delta));
    print(paste('Volume(SG) = ',vol1));
    print(paste('Volume(C+1)= ',vol2));
    print(paste(c1));
    print(out);
    #print(paste(c(num,out)));
    print(paste(c2));
    #print(paste(c(num,sg)));
    print(sg);
  }
  return(salida);
}