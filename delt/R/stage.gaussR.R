stage.gaussR<-function(dendat,M,mugrid,siggrid=1,sigeka=TRUE,sampstart=TRUE)
{
n<-length(dendat)
dict.card<-length(mugrid)
dict.card.sig<-length(siggrid)

muut<-matrix(0,M,1)      #gives the means of the final mixture
sigit<-matrix(0,M,1)     #gives the std:s of the final mixture
piit<-matrix(0,M-1,1)
for (i in 1:(M-1)) piit[i]<-2/(i+2)
riskit<-matrix(0,dict.card,dict.card.sig)

if (sampstart){
   muut[1]<-mean(dendat)
   sigit[1]<-sqrt(var(dendat))
}
else{
   # haetaan 1. termi
   for (i in 1:dict.card){
      for (ii in 1:dict.card.sig){
         sqint<-gaussprod(0,0,siggrid[ii],siggrid[ii])
         val<-0
         for (j in 1:n){   
             point<-dendat[j]
             evapoint<-(point-mugrid[i])/siggrid[ii]
             val<-val+evanor(evapoint)/siggrid[ii]
         }
         riskit[i,ii]<--2*val+sqint
      }
   }
   mind<-which.min(t(riskit))
   sarat<-dict.card.sig
   imin<-ceiling(mind/sarat)
   jmin<-mind-(imin-1)*sarat

   muut[1]<-mugrid[imin]
   if (sigeka) sigit[1]<-1 else sigit[1]<-siggrid[jmin]
}

# haetaan termit 2-M
curmix<-matrix(0,M,1)   #estimaatin painot
curmix[1]<-1            #alussa vain yksi simppeli funktio
k<-1
while (k <= (M-1)){
   for (i in 1:dict.card){
      for (ii in 1:dict.card.sig){
          # calculate the -2*average of evaluations
          val<-0
          for (j in 1:n){   
              point<-dendat[j]
              evapoint<-(point-mugrid[i])/siggrid[ii]
              val<-val+evanor(evapoint)/siggrid[ii]
          }
          # calculate the inner product of the candidate with the k-1 estimate
          prodint<-0
          jj<-1
          while (jj<=k){
             prodint<-prodint+curmix[jj]*gaussprod(muut[jj],mugrid[i],
                                         sigit[jj],siggrid[ii])
             jj<-jj+1
          }
          # calculate the risk at the k:th step
          gammanpik<--2*piit[k]*val/n+piit[k]^2*gaussprod(0,0,siggrid[ii],
                                                siggrid[ii])
          riskit[i,ii]<-gammanpik+2*(1-piit[k])*piit[k]*prodint
      }
   }  
   mind<-which.min(t(riskit))
   sarat<-dict.card.sig
   imin<-ceiling(mind/sarat)
   jmin<-mind-(imin-1)*sarat

   muut[k+1]<-mugrid[imin]
   sigit[k+1]<-siggrid[jmin]

   curmix[1:k]<-(1-piit[k])*curmix[1:k]
   curmix[k+1]<-piit[k]
   k<-k+1
}

#sig<-matrix(1,M,1)
#et<-eval.func("mixt",N,sig=sigit,M=muut,p=curmix)  

return(list(muut=muut,sigit=sigit,curmix=curmix))
}  


                            

