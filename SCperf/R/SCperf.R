SCperf <-
function(phi,theta,L=L,SL=0.95)
{   if (L==0){cat("L is at least the review period which is one, ...\n")
             }
    else 
  {
   #Calculating the variance of the demand using ARMAtoMA

   arma<-ARMAtoMA(ar=phi, ma=theta, 1000);
   VarD<-sum(arma^2)+1;
              
    #Calculating the bullwhip effect

     values<-ARMAtoMA(ar=phi, ma=theta, L);

       total = 0      #Calculating de duble sum in the formula of BE
       for (i in 1:L)
       {    valsum<-  sum (values[i:L]);
            if( i==1)  {total <- valsum;}
            else  {total <- total + values[i-1] * valsum;}
       }
       be<-1+2*total/VarD;
       
      t=ifelse(phi==0,1,be)
            
     #Calculating the variance during the LT
      arma1<-ARMAtoMA(ar=phi, ma=theta,L);
      arma2<-c(1,arma1);

      totalLT <- 0
      for (i in 1:L)
       { valsumLT<-  (sum(arma2[1:i]))^2;
         totalLT<- totalLT + valsumLT;
       }

        VarLT<-totalLT;

        #Calculating de SS=z*sigma*L^0.5 and SSLT=z*sigmaLT
        z<- qnorm(p=SL, mean = 0, sd = 1)
        SS<-z*(VarD*L)^0.5
        SSLT<-z*(VarLT^0.5)
   }
     f<-c(M=t[1],VarD=VarD,VarDL=VarLT,SS=SS,SSL=SSLT,z=z)
     return(f)
}

