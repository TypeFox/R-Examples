.moments<-function(c,lambda){
   a=lambda+c;
   b=lambda-c;
   if(b<0){ 
      b=0;
   }
   if(b>0){ 
      # den=poisson(lambda,a)-poisson(lambda,b-1);
      den=ppois(a,lambda)-ppois(b-1,lambda);
   }
   if(b==0){
      # den=poisson(lambda,a);
      den=ppois(a,lambda);
   }
   mu=mat.or.vec(4,1); 
   # mom es global y se usa fuera de esta función
   mom=mat.or.vec(5,1); 
   for(r in 1:4){
      poisA=0;
      poisB=0;
      if((a-r) >=0){ poisA=ppois(a,lambda)-ppois(a-r,lambda); }
      if((a-r) < 0){ poisA=ppois(a,lambda); }
      if((b-r-1) >=0){ poisB=ppois(b-1,lambda)-ppois(b-r-1,lambda); }
      if((b-r-1) < 0 && (b-1)>=0){ poisB=ppois(b-1,lambda); }
      if((b-r-1) < 0 && (b-1) < 0){ poisB=0; }
      mu[r]=(lambda^r)*(1-(poisA-poisB)/den);
   }
   mom[1]=mu[1];
   mom[2]=mu[2]+mu[1]-mu[1]^2;
   mom[3]=mu[3]+mu[2]*(3-3*mu[1])+(mu[1]-3*mu[1]^2+2*mu[1]^3);
   mom[4]=mu[4]+mu[3]*(6-4*mu[1])+mu[2]*(7-12*mu[1]+6*mu[1]^2)+mu[1]-4*mu[1]^2+6*mu[1]^3-3*mu[1]^4;
   mom[5]=den;
   mom
}
.truncpoi<-function(c,x,n,k){
   m=matrix(0,k,5);
   for(i in 1:k){
    lambda=x[i];
    mom = .moments(c,lambda);
    for(j in 1:5){
     m[i,j]=mom[j];
    }
   }
   for(i in 1:k){
    m[i,4]=m[i,4]-3*m[i,2]^2;
   }
   
   #s1=m[+,1];
   #s2=m[+,2];
   #s3=m[+,3];
   #s4=m[+,4];
   s=colSums(m);
   s1=s[1];
   s2=s[2];
   s3=s[3];
   s4=s[4];

   probn=1/(ppois(n,n)-ppois(n-1,n));
   z=(n-s1)/sqrt(s2);
   g1=s3/(s2^(3/2));
   g2=s4/(s2^2);
   poly=1+g1*(z^3-3*z)/6+g2*(z^4-6*z^2+3)/24
         +g1^2*(z^6-15*z^4+45*z^2-15)/72;
   f=poly*exp(-z^2/2)/(sqrt(2)*gamma(0.5));
   probx=1;
   for(i in 1:k){
    probx=probx*m[i,5];
   }
   return(probn*probx*f/sqrt(s2));
}